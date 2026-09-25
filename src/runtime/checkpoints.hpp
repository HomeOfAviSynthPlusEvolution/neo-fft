#pragma once
#include "kernels/kalman.hpp"
#include <array>
#include <memory>
#include <mutex>
namespace neo_fft::runtime {
struct Checkpoint {
  int frame = 0;
  std::size_t control_bytes = 0;
  std::array<KalmanState,4> planes;
  std::size_t bytes() const {
    std::size_t n=add_size(sizeof(*this),control_bytes);
    for (const auto& p:planes) n=add_size(n,p.bytes()-sizeof(p));
    return n;
  }
};
// std::shared_ptr rebinds this allocator to its actual control-block type.
// Only allocation touches the counter; deallocation follows object destruction.
template<class T> struct ControlAllocator {
  using value_type=T;
  std::size_t* bytes;
  explicit ControlAllocator(std::size_t* b):bytes(b) {}
  template<class U> ControlAllocator(const ControlAllocator<U>& a):bytes(a.bytes) {}
  T* allocate(std::size_t count) {
    T* p=std::allocator<T>{}.allocate(count);*bytes=mul_size(count,sizeof(T));return p;
  }
  void deallocate(T* p,std::size_t count) {std::allocator<T>{}.deallocate(p,count);}
  template<class U> bool operator==(const ControlAllocator<U>& a) const {return bytes==a.bytes;}
  template<class U> bool operator!=(const ControlAllocator<U>& a) const {return !(*this==a);}
};
// Fixed metadata, no clip-sized map. Leases keep evicted immutable states alive.
class Checkpoints {
public:
  using Lease=std::shared_ptr<const Checkpoint>;
  Checkpoints():budget_(64*1024*1024),retain_one_(true) {}
  explicit Checkpoints(std::size_t budget):budget_(budget),retain_one_(false) {}
  Lease acquire(int frame,int earliest=0) {
    std::lock_guard<std::mutex> lock(mutex_);
    Entry* best=nullptr;
    for (auto& e:entries_) if(e.value && e.value->frame>=earliest && e.value->frame<=frame &&
                            (!best || e.value->frame>best->value->frame)) best=&e;
    if(!best) return {};
    touch(*best); return best->value;
  }
  void publish(std::unique_ptr<Checkpoint> owned) {
    if(!owned || owned->frame<1) return;
    auto* raw=owned.release();
    Lease value(raw,std::default_delete<Checkpoint>{},ControlAllocator<Checkpoint>(&raw->control_bytes));
    const auto bytes=value->bytes();
    // Large geometries must retain a complete state to continue sequentially
    // instead of restarting warmup on every request. Explicit budgets are hard.
    const auto required=add_size(sizeof(*this),bytes);
    const auto limit=retain_one_ && required>budget_ ? required : budget_;
    if(required>limit) return;
    std::lock_guard<std::mutex> lock(mutex_);
    // Histories can differ after a seek: keep the first published snapshot
    // while retained. Never replace a state already leased by another request.
    for(auto& e:entries_) if(e.value && e.value->frame==value->frame) {touch(e);return;}
    while(used_>limit-required || count()==entries_.size()) {
      Entry* oldest=nullptr;
      for(auto& e:entries_) if(e.value && (!oldest || e.age>oldest->age)) oldest=&e;
      used_-=oldest->value->bytes(); oldest->value.reset();
    }
    for(auto& e:entries_) if(!e.value) { e.value=std::move(value); e.age=4; used_+=bytes;touch(e);break; }
  }
  std::size_t bytes() const {std::lock_guard<std::mutex> lock(mutex_);return used_+sizeof(*this);}
private:
  struct Entry {Lease value;unsigned age=0;};
  void touch(Entry& selected) {
    for(auto& e:entries_) if(e.value && &e!=&selected && e.age<selected.age) ++e.age;
    selected.age=0;
  }
  std::size_t count() const {std::size_t n=0;for(const auto& e:entries_) n+=bool(e.value);return n;}
  std::array<Entry,4> entries_{};
  const std::size_t budget_;
  const bool retain_one_;
  std::size_t used_=0;
  mutable std::mutex mutex_;
};
} // namespace neo_fft::runtime
