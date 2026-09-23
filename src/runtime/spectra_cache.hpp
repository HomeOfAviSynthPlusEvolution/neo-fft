#pragma once
#include "runtime/checkpoints.hpp"
#include <condition_variable>
#include <exception>
#include <list>

namespace neo_fft::runtime {
// Instance-local, immutable raw spatial spectra. One entry reserves all selected
// planes of one frame, including planes not built yet. Pinned entries cannot be
// evicted: reservations therefore include builders and every outstanding lease.
class SpectraCache {
  struct Frame {
    explicit Frame(int n):number(n) {}
    int number;
    std::size_t charge=0;
    std::array<std::unique_ptr<std::complex<float>[]>,3> data;
    enum Status { Empty, Building, Ready, Failed };
    std::array<Status,3> status{};
    std::array<std::exception_ptr,3> errors;
    std::condition_variable changed;
  };
  using FramePtr=std::shared_ptr<Frame>;
  using Entries=std::list<FramePtr,ControlAllocator<FramePtr>>;
public:
  using Lease=std::shared_ptr<const std::complex<float>>;
  SpectraCache(std::array<std::size_t,3> bins,int temporal,int frames=-1,int mb=512)
      :bins_(bins),temporal_(temporal),requested_frames_(frames),
       budget_(mul_size(std::size_t(mb==-1 ? 512 : std::max(0,mb)),1024*1024)),
       entries_(ControlAllocator<FramePtr>(&node_bytes_)) {
    require(frames>=-1 && mb>=-1,"FFT3D cache limits must be >= -1");
    require(temporal>=1 && temporal<=5,"invalid cache temporal size");
    for(auto n:bins_)payload_=add_size(payload_,mul_size(n,sizeof(std::complex<float>)));
  }
  void inform_threads(int threads) {
    if(threads<=0)return;
    std::lock_guard<std::mutex> lock(mutex_);threads_=threads;
    trim();
  }
  std::size_t frame_limit() const {
    std::lock_guard<std::mutex> lock(mutex_);return limit();
  }
  std::size_t bytes() const {std::lock_guard<std::mutex> lock(mutex_);return used_;}
  std::size_t peak_bytes() const {std::lock_guard<std::mutex> lock(mutex_);return peak_;}
  std::size_t builds() const {std::lock_guard<std::mutex> lock(mutex_);return builds_;}
  std::size_t hits() const {std::lock_guard<std::mutex> lock(mutex_);return hits_;}
  std::size_t waits() const {std::lock_guard<std::mutex> lock(mutex_);return waits_;}

  // The builder runs synchronously, performs no host fetch or nested cache
  // acquisition, and publishes a complete finite plane. Thus same-key waiting
  // cannot form request cycles. Budget exhaustion never waits for leases.
  template<class Build> Lease get(int frame,int plane,Build&& build) {
    require(frame>=0 && plane>=0 && plane<3,"invalid spectrum cache key");
    std::unique_lock<std::mutex> lock(mutex_);
    trim();
    if(!budget_ || !limit() || temporal_<=1 || !bins_[plane])return {};
    FramePtr entry;
    for(auto it=entries_.begin();it!=entries_.end();++it)if((*it)->number==frame) {
      entry=*it;entries_.splice(entries_.end(),entries_,it);break;
    }
    if(!entry) {
      // Measure actual shared ownership and list-node allocations. Only this
      // small candidate metadata precedes admission; no spectra are allocated.
      std::size_t control=0;
      entry=std::allocate_shared<Frame>(ControlAllocator<Frame>(&control),frame);
      entries_.push_back(entry);
      try {entry->charge=add_size(payload_,add_size(control,node_bytes_));}
      catch(...) {entries_.pop_back();throw;}
      if(entry->charge>budget_) {entries_.pop_back();return {};}
      while(used_>budget_-entry->charge || entries_.size()>limit()) {
        if(!evict_one()) {entries_.pop_back();return {};}
      }
      used_+=entry->charge;peak_=std::max(peak_,used_);
    }
    if(entry->status[plane]==Frame::Building) {
      ++waits_;
      entry->changed.wait(lock,[&]{return entry->status[plane]!=Frame::Building;});
      if(entry->status[plane]==Frame::Failed)std::rethrow_exception(entry->errors[plane]);
    }
    if(entry->status[plane]==Frame::Ready) {++hits_;return Lease(entry,entry->data[plane].get());}
    // A failed producer never poisons future requests: use the uncached path
    // until this entry is evicted. Existing waiters receive its original error.
    if(entry->status[plane]==Frame::Failed)return {};
    entry->status[plane]=Frame::Building;++builds_;
    lock.unlock();
    try {
      auto data=std::make_unique<std::complex<float>[]>(bins_[plane]);
      build(data.get());
      lock.lock();entry->data[plane]=std::move(data);entry->status[plane]=Frame::Ready;
      entry->changed.notify_all();return Lease(entry,entry->data[plane].get());
    } catch(...) {
      if(!lock.owns_lock())lock.lock();
      entry->errors[plane]=std::current_exception();entry->status[plane]=Frame::Failed;
      entry->changed.notify_all();throw;
    }
  }
private:
  std::size_t limit() const {
    return requested_frames_==-1 ? std::size_t(temporal_)+std::size_t(threads_)-1 : std::size_t(requested_frames_);
  }
  bool evict_one() {
    for(auto it=entries_.begin();it!=entries_.end();++it)if(it->use_count()==1) {
      used_-=(*it)->charge;entries_.erase(it);return true;
    }
    return false;
  }
  void trim() {while(entries_.size()>limit() && evict_one()) {}}
  const std::array<std::size_t,3> bins_;
  const int temporal_,requested_frames_;
  const std::size_t budget_;
  std::size_t payload_=0,node_bytes_=0,used_=0,peak_=0,builds_=0,hits_=0,waits_=0;
  int threads_=1;
  mutable std::mutex mutex_;
  Entries entries_;
};
} // namespace neo_fft::runtime
