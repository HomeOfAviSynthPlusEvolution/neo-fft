#pragma once
#include "base/checked.hpp"
#include <memory>
#include <mutex>

namespace neo_fft::runtime {
// Build privately without this lock (in particular, never fetch host frames under
// it). Concurrent duplicate builds are allowed; only complete immutable models
// can be published. Throwing builders leave the cache absent and retryable.
class PublishedModel {
public:
  using Model = std::shared_ptr<const std::vector<float>>;
  Model get() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return model_;
  }
  Model publish(std::vector<float> values) {
    for (float v : values) require(std::isfinite(v) && v >= 0, "invalid sampled power");
    auto complete = std::make_shared<const std::vector<float>>(std::move(values));
    return publish(std::move(complete));
  }
  Model publish(Model complete) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!model_) model_ = std::move(complete);
    return model_;
  }
private:
  mutable std::mutex mutex_;
  Model model_;
};
} // namespace neo_fft::runtime
