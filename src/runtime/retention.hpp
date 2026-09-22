#pragma once
#include <cstddef>
#include <mutex>
namespace neo_fft::runtime {
// One admission counter shared by every plane's idle workspace pool.
class Retention {
public:
  explicit Retention(std::size_t count,std::size_t bytes=64*1024*1024):max_count_(count),max_bytes_(bytes) {}
  bool acquire(std::size_t bytes) {
    std::lock_guard<std::mutex> lock(mutex_);
    if(count_==max_count_ || bytes>max_bytes_-bytes_) return false;
    ++count_;bytes_+=bytes;return true;
  }
  void release(std::size_t bytes) {
    std::lock_guard<std::mutex> lock(mutex_);--count_;bytes_-=bytes;
  }
  std::size_t bytes() const {std::lock_guard<std::mutex> lock(mutex_);return bytes_;}
  std::size_t count() const {std::lock_guard<std::mutex> lock(mutex_);return count_;}
private:
  const std::size_t max_count_,max_bytes_;
  std::size_t count_=0,bytes_=0;
  mutable std::mutex mutex_;
};
} // namespace neo_fft::runtime
