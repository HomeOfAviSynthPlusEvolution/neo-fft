#include "runtime/workspace_pool.hpp"
#include "base/checked.hpp"

namespace neo_fft::runtime {

WorkspacePool::WorkspacePool(WorkspaceBudget budget, std::size_t max_capacity, std::shared_ptr<Retention> retention)
    : budget_(budget), max_capacity_(max_capacity), retention_(retention ? std::move(retention) : std::make_shared<Retention>(max_capacity)) {
  require(max_capacity > 0, "max_capacity must be positive");
  idle_.reserve(max_capacity);
}

WorkspacePool::~WorkspacePool() {
  for(const auto& ws:idle_) retention_->release(ws->retained_bytes());
}

WorkspaceLease WorkspacePool::acquire() {
  std::unique_lock<std::mutex> lock(mutex_);
  if (!idle_.empty()) {
    auto ws = std::move(idle_.back());
    idle_.pop_back();
    retention_->release(ws->retained_bytes());
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
    throw;
  }
}

void WorkspacePool::release(std::unique_ptr<Workspace> ws) noexcept {
  if (!ws) {
    return;
  }
  std::lock_guard<std::mutex> lock(mutex_);
  if (idle_.size()<max_capacity_ && retention_->acquire(ws->retained_bytes())) idle_.push_back(std::move(ws));
  --active_count_;
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
