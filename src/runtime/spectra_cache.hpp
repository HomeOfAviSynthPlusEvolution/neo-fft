#pragma once
#include "runtime/checkpoints.hpp"
#include <condition_variable>
#include <exception>
#include <map>
#include <tuple>

namespace neo_fft::runtime {
// Immutable spatial block rows. Builders and consumers pin the same entry;
// every live payload remains charged until the last lease permits eviction.
class SpectraCache {
  using Key=std::tuple<int,int,int>;
  struct Row {
    explicit Row(Key k):key(k) {}
    Key key;
    std::size_t charge=0;
    std::unique_ptr<std::complex<float>[]> data;
    enum Status { Building, Ready, Failed } status=Building;
    std::exception_ptr error;
    std::condition_variable changed;
    Row* previous=nullptr;
    Row* next=nullptr;
  };
  using RowPtr=std::shared_ptr<Row>;
  using Item=std::pair<const Key,RowPtr>;
  using Index=std::map<Key,RowPtr,std::less<Key>,ControlAllocator<Item>>;
  struct Frame {std::size_t rows=0;std::uint64_t age=0;};
  using FrameItem=std::pair<const int,Frame>;
  using Frames=std::map<int,Frame,std::less<int>,ControlAllocator<FrameItem>>;
public:
  using Lease=std::shared_ptr<const std::complex<float>>;
  struct Stats {
    std::uint64_t requests=0,hits=0,waits=0,builds=0,bypasses=0,evictions=0;
    std::uint64_t requested_bins=0,computed_bins=0;
    std::size_t bytes=0,peak_bytes=0,rows=0,frames=0;
  };
  SpectraCache(std::array<std::size_t,3> bins,int temporal,int frames=-1,int mb=512,
               std::array<int,3> rows={1,1,1})
      :bins_(bins),rows_(rows),temporal_(temporal),requested_frames_(frames),
       budget_(mul_size(std::size_t(mb==-1 ? 512 : std::max(0,mb)),1024*1024)),
       index_(ControlAllocator<Item>(&node_bytes_)),frames_(ControlAllocator<FrameItem>(&frame_bytes_)) {
    require(frames>=-1 && mb>=-1,"FFT3D cache limits must be >= -1");
    require(temporal>=1 && temporal<=5,"invalid cache temporal size");
    for(int p=0;p<3;++p) {
      require(rows_[p]>=0 && (!bins_[p] || rows_[p]>0),"invalid spectrum cache shape");
      payload_[p]=mul_size(bins_[p],sizeof(std::complex<float>));
    }
  }
  void inform_threads(int threads) {
    if(threads<=0)return;
    std::lock_guard<std::mutex> lock(mutex_);threads_=threads;trim();
  }
  std::size_t frame_limit() const {std::lock_guard<std::mutex> lock(mutex_);return limit();}
  Stats stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto result=stats_;result.bytes=used_;result.rows=index_.size();result.frames=frames_.size();return result;
  }
  std::size_t bytes() const {return stats().bytes;}
  std::size_t peak_bytes() const {return stats().peak_bytes;}
  std::uint64_t builds() const {return stats().builds;}
  std::uint64_t hits() const {return stats().hits;}
  std::uint64_t waits() const {return stats().waits;}

  // Synchronous builders must not fetch host frames or acquire another cache
  // entry. Thus waiting for the same row cannot form dependency cycles.
  // Exhausted budgets fall back immediately, without waiting for lease release.
  template<class Build> Lease get(int frame,int plane,int row,Build&& build) {
    require(frame>=0 && plane>=0 && plane<3 && row>=0 && row<rows_[plane],"invalid spectrum cache key");
    std::unique_lock<std::mutex> lock(mutex_);
    trim();++stats_.requests;stats_.requested_bins+=bins_[plane];
    const auto bypass=[&]() -> Lease {++stats_.bypasses;stats_.computed_bins+=bins_[plane];return {};};
    if(!budget_ || !limit() || temporal_<=1 || !bins_[plane] || payload_[plane]>budget_)return bypass();
    const Key key{frame,plane,row};
    auto found=index_.find(key);
    if(found!=index_.end()) {
      auto entry=found->second;touch(entry.get());frames_.find(frame)->second.age=++age_;
      if(entry->status==Row::Building) {
        ++stats_.waits;
        entry->changed.wait(lock,[&]{return entry->status!=Row::Building;});
        if(entry->status==Row::Failed)std::rethrow_exception(entry->error);
      }
      if(entry->status==Row::Failed)return bypass();
      ++stats_.hits;return Lease(entry,entry->data.get());
    }
    // Only small metadata precedes admission. Measure the implementation's
    // actual shared-control and map-node allocation sizes, including frame keys.
    std::size_t control=0;
    auto entry=std::allocate_shared<Row>(ControlAllocator<Row>(&control),key);
    auto frame_entry=frames_.try_emplace(frame);
    try {
      index_.emplace(key,entry);
      entry->charge=add_size(payload_[plane],add_size(control,node_bytes_));
      const auto charge=add_size(entry->charge,frame_entry.second ? frame_bytes_ : 0);
      if(charge>budget_ || charge>std::numeric_limits<std::size_t>::max()-used_) {
        index_.erase(key);if(frame_entry.second)frames_.erase(frame_entry.first);
        return bypass();
      }
      used_+=charge;
    } catch(...) {
      index_.erase(key);if(frame_entry.second)frames_.erase(frame_entry.first);throw;
    }
    ++frame_entry.first->second.rows;frame_entry.first->second.age=++age_;touch(entry.get());
    while(frames_.size()>limit())if(!evict_frame()) {erase(entry.get());return bypass();}
    while(used_>budget_)if(!evict_row()) {erase(entry.get());return bypass();}
    stats_.peak_bytes=std::max(stats_.peak_bytes,used_);
    ++stats_.builds;stats_.computed_bins+=bins_[plane];
    lock.unlock();
    try {
      auto data=std::make_unique<std::complex<float>[]>(bins_[plane]);build(data.get());
      lock.lock();entry->data=std::move(data);entry->status=Row::Ready;
      entry->changed.notify_all();return Lease(entry,entry->data.get());
    } catch(...) {
      if(!lock.owns_lock())lock.lock();
      entry->error=std::current_exception();entry->status=Row::Failed;
      entry->changed.notify_all();throw;
    }
  }
private:
  std::size_t limit() const {
    return requested_frames_==-1 ? std::size_t(temporal_)+std::size_t(threads_)-1 : std::size_t(requested_frames_);
  }
  void unlink(Row* entry) {
    if(entry->previous)entry->previous->next=entry->next;else if(first_==entry)first_=entry->next;
    if(entry->next)entry->next->previous=entry->previous;else if(last_==entry)last_=entry->previous;
    entry->previous=entry->next=nullptr;
  }
  void touch(Row* entry) {
    unlink(entry);entry->previous=last_;
    if(last_)last_->next=entry;else first_=entry;
    last_=entry;
  }
  void erase(Row* entry) {
    const auto key=entry->key;unlink(entry);used_-=entry->charge;
    auto frame=frames_.find(std::get<0>(key));
    if(--frame->second.rows==0) {used_-=frame_bytes_;frames_.erase(frame);}
    index_.erase(key);
  }
  bool evict_row() {
    for(auto* row=first_;row;row=row->next)if(index_.find(row->key)->second.use_count()==1) {
      erase(row);++stats_.evictions;return true;
    }
    return false;
  }
  bool evict_frame() {
    auto victim=frames_.end();
    for(auto frame=frames_.begin();frame!=frames_.end();++frame) {
      if(victim!=frames_.end() && frame->second.age>=victim->second.age)continue;
      bool pinned=false;
      for(auto it=index_.lower_bound(Key{frame->first,0,0});it!=index_.end() && std::get<0>(it->first)==frame->first;++it)
        if(it->second.use_count()!=1) {pinned=true;break;}
      if(!pinned)victim=frame;
    }
    if(victim==frames_.end())return false;
    const int number=victim->first;
    auto it=index_.lower_bound(Key{number,0,0});
    while(it!=index_.end() && std::get<0>(it->first)==number) {
      auto* row=(it++)->second.get();erase(row);++stats_.evictions;
    }
    return true;
  }
  void trim() {while(frames_.size()>limit() && evict_frame()) {}}
  const std::array<std::size_t,3> bins_;
  const std::array<int,3> rows_;
  const int temporal_,requested_frames_;
  const std::size_t budget_;
  std::array<std::size_t,3> payload_{};
  std::size_t node_bytes_=0,frame_bytes_=0,used_=0;
  std::uint64_t age_=0;
  Stats stats_;
  int threads_=1;
  mutable std::mutex mutex_;
  Index index_;
  Frames frames_;
  Row* first_=nullptr;
  Row* last_=nullptr;
};
} // namespace neo_fft::runtime
