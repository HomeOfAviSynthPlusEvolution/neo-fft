#pragma once
#include "base/checked.hpp"
#include <atomic>
#include <future>
#include <memory>
namespace neo_fft::runtime {
// Shared per instance; busy instances fall back to their own calling thread.
// Each invocation joins only workers it created, including during unwinding.
class Executor {
public:
  explicit Executor(int workers=1):workers_(workers),available_(workers-1) {require(workers>0 && workers<=16,"invalid worker capacity");}
  int workers() const noexcept {return workers_;}
  int peak() const noexcept {return peak_.load();}
  template<class Function> void run(int count,Function&& function) const {
    if(count<=0) return;
    std::atomic<int> next{0};
    auto consume=[&] {for(;;) {
      const int index=next.fetch_add(1);if(index>=count) break;
      const int active=active_.fetch_add(1)+1;int peak=peak_.load();
      while(peak<active && !peak_.compare_exchange_weak(peak,active)) {}
      struct Active {std::atomic<int>& value;~Active(){value.fetch_sub(1);}} guard{active_};
      function(index);
    }};
    std::vector<std::future<void>> jobs;jobs.reserve(std::size_t(std::min(count-1,workers_-1)));
    std::exception_ptr failure;
    for(int i=1;i<count && i<workers_;++i) {
      int available=available_.load();
      while(available>0 && !available_.compare_exchange_weak(available,available-1)) {}
      if(available==0) break;
      try {
        jobs.push_back(std::async(std::launch::async,[&] {
          struct Slot {const Executor* executor;~Slot(){executor->available_.fetch_add(1);}} slot{this};
          consume();
        }));
      } catch(...) {available_.fetch_add(1);failure=std::current_exception();break;}
    }
    try {consume();} catch(...) {failure=std::current_exception();}
    for(auto& job:jobs) try {job.get();} catch(...) {if(!failure)failure=std::current_exception();}
    if(failure) std::rethrow_exception(failure);
  }
private:
  int workers_;
  mutable std::atomic<int> available_,active_{0},peak_{0};
};
} // namespace neo_fft::runtime
