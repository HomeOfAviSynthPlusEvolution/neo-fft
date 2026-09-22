#include "runtime/workspace_pool.hpp"
#include "base/checked.hpp"

namespace neo_fft::runtime {

WorkspacePool::WorkspacePool(WorkspaceBudget budget, std::size_t max_capacity)
    : budget_(budget), max_capacity_(max_capacity) {
  require(max_capacity > 0, "max_capacity must be positive");
  idle_.reserve(max_capacity);
}

WorkspaceLease WorkspacePool::acquire() {
  std::unique_lock<std::mutex> lock(mutex_);
  cv_.wait(lock, [this] { return !idle_.empty() || active_count_ < max_capacity_; });
  if (!idle_.empty()) {
    auto ws = std::move(idle_.back());
    idle_.pop_back();
    ++active_count_;
    return WorkspaceLease(this, std::move(ws));
  }
  ++active_count_;
  lock.unlock();
  try {
    auto ws = std::make_unique<Workspace>(budget_);
    return WorkspaceLease(this, std::move(ws));
  } catch (...) {
    lock.lock();
    --active_count_;
    cv_.notify_one();
    throw;
  }
}

void WorkspacePool::release(std::unique_ptr<Workspace> ws) noexcept {
  if (!ws) {
    return;
  }
  std::lock_guard<std::mutex> lock(mutex_);
  idle_.push_back(std::move(ws));
  --active_count_;
  cv_.notify_one();
}

std::size_t WorkspacePool::active_count() const noexcept {
  std::lock_guard<std::mutex> lock(mutex_);
  return active_count_;
}

std::size_t WorkspacePool::idle_count() const noexcept {
  std::lock_guard<std::mutex> lock(mutex_);
  return idle_.size();
}

WorkspaceLease::~WorkspaceLease() {
  reset();
}

void WorkspaceLease::reset() noexcept {
  if (pool_ && ws_) {
    pool_->release(std::move(ws_));
  }
  pool_ = nullptr;
  ws_.reset();
}

} // namespace neo_fft::runtime
