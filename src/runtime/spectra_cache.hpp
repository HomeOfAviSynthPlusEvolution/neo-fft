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
  static constexpr int default_mb=128;
  // Request storage belongs to its host request, not to cached row entries.
  // The cache must outlive registrations; no host calls occur under its mutex.
  class Request {
    friend class SpectraCache;
    Request(SpectraCache& cache,int first,int last,std::array<bool,3> planes)
        :cache_(cache),first_(first),last_(last) {
      for(int p=0;p<3;++p)next_row_[p]=planes[p] ? 0 : cache.rows_[p];
    }
  public:
    Request(const Request&)=delete;
    Request& operator=(const Request&)=delete;
    ~Request() {
      if(!linked_)return;
      std::lock_guard<std::mutex> lock(cache_.mutex_);
      if(previous_)previous_->next_=next_;else cache_.requests_=next_;
      if(next_)next_->previous_=previous_;
      --cache_.stats_.active_requests;
    }
    void complete_row(int plane,int row) {
      std::lock_guard<std::mutex> lock(cache_.mutex_);
      next_row_[plane]=row+1;
    }
  private:
    SpectraCache& cache_;
    int first_,last_;
    std::array<int,3> next_row_{};
    bool linked_=false;
    Request* previous_=nullptr;
    Request* next_=nullptr;
  };
  using Lease=std::shared_ptr<const std::complex<float>>;
  struct Stats {
    std::uint64_t requests=0,hits=0,waits=0,builds=0,bypasses=0,evictions=0;
    std::uint64_t requested_bins=0,computed_bins=0;
    std::size_t bytes=0,peak_bytes=0,rows=0,frames=0;
    std::uint64_t registrations=0,active_requests=0,peak_requests=0;
    std::uint64_t evicted_needed_rows=0,evicted_consumers=0;
  };
  SpectraCache(std::array<std::size_t,3> bins,int temporal,int frames=-1,int mb=default_mb,
               std::array<int,3> rows={1,1,1})
      :bins_(bins),rows_(rows),temporal_(temporal),requested_frames_(frames),
       budget_(mul_size(std::size_t(mb==-1 ? default_mb : std::max(0,mb)),1024*1024)),
       index_(ControlAllocator<Item>(&node_bytes_)),frames_(ControlAllocator<FrameItem>(&frame_bytes_)) {
    require(frames>=-1 && mb>=-1,"FFT3D cache limits must be >= -1");
    require(temporal>=1 && temporal<=5,"invalid cache temporal size");
    for(int p=0;p<3;++p) {
      require(rows_[p]>=0 && (!bins_[p] || rows_[p]>0),"invalid spectrum cache shape");
      payload_[p]=mul_size(bins_[p],sizeof(std::complex<float>));
    }
  }
  std::unique_ptr<Request> register_request(int first,int last,std::array<bool,3> planes={true,true,true}) {
    require(first>=0 && last>=first && last-first<temporal_,"invalid cache request interval");
    if(!budget_ || requested_frames_==0 || temporal_<=1)return {};
    auto request=std::unique_ptr<Request>(new Request(*this,first,last,planes));
    std::lock_guard<std::mutex> lock(mutex_);
    request->next_=requests_;if(requests_)requests_->previous_=request.get();requests_=request.get();
    request->linked_=true;registered_=true;++stats_.registrations;++stats_.active_requests;
    stats_.peak_requests=std::max(stats_.peak_requests,stats_.active_requests);
    return request;
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
  std::uint64_t pending_consumers(const Key& key) const {
    const auto [frame,plane,row]=key;
    std::uint64_t count=0;
    for(auto* request=requests_;request;request=request->next_)
      if(frame>=request->first_ && frame<=request->last_ && row>=request->next_row_[plane])++count;
    return count;
  }
  void evict(Row* row) {
    const auto pending=pending_consumers(row->key);
    stats_.evicted_needed_rows+=pending!=0;stats_.evicted_consumers+=pending;
    erase(row);++stats_.evictions;
  }
  bool evict_row() {
    if(registered_) {
      // The ordered index prefers earlier source frames within each demand
      // class. Pending consumers are a soft priority, never an extra pin.
      Row* needed=nullptr;
      for(const auto& item:index_)if(item.second.use_count()==1) {
        auto* row=item.second.get();
        if(!pending_consumers(row->key)) {evict(row);return true;}
        if(!needed)needed=row;
      }
      if(needed) {evict(needed);return true;}
    } else {
      for(auto* row=first_;row;row=row->next)if(index_.find(row->key)->second.use_count()==1) {
        evict(row);return true;
      }
    }
    return false;
  }
  bool evict_frame() {
    auto victim=frames_.end();bool victim_needed=true;
    for(auto frame=frames_.begin();frame!=frames_.end();++frame) {
      if(!registered_ && victim!=frames_.end() && frame->second.age>=victim->second.age)continue;
      bool pinned=false,needed=false;
      for(auto it=index_.lower_bound(Key{frame->first,0,0});it!=index_.end() && std::get<0>(it->first)==frame->first;++it) {
        if(it->second.use_count()!=1) {pinned=true;break;}
        if(registered_ && pending_consumers(it->first))needed=true;
      }
      if(pinned)continue;
      if(!registered_ || victim==frames_.end() || (victim_needed && !needed)) {
        victim=frame;victim_needed=needed;
      }
      if(registered_ && !needed)break;
    }
    if(victim==frames_.end())return false;
    const int number=victim->first;
    auto it=index_.lower_bound(Key{number,0,0});
    while(it!=index_.end() && std::get<0>(it->first)==number) {
      auto* row=(it++)->second.get();evict(row);
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
  Request* requests_=nullptr;
  bool registered_=false;
  Row* first_=nullptr;
  Row* last_=nullptr;
};
} // namespace neo_fft::runtime
