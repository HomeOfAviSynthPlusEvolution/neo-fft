#pragma once
#include "runtime/workspace.hpp"
#include <condition_variable>
#include <memory>
#include <mutex>
#include <vector>

namespace neo_fft::runtime {

class WorkspacePool;

class WorkspaceLease {
public:
  WorkspaceLease() noexcept = default;
  WorkspaceLease(WorkspacePool* pool, std::unique_ptr<Workspace> ws) noexcept
      : pool_(pool), ws_(std::move(ws)) {}
  ~WorkspaceLease();

  WorkspaceLease(WorkspaceLease&& other) noexcept
      : pool_(other.pool_), ws_(std::move(other.ws_)) {
    other.pool_ = nullptr;
  }

  WorkspaceLease& operator=(WorkspaceLease&& other) noexcept {
    if (this != &other) {
      reset();
      pool_ = other.pool_;
      ws_ = std::move(other.ws_);
      other.pool_ = nullptr;
    }
    return *this;
  }

  WorkspaceLease(const WorkspaceLease&) = delete;
  WorkspaceLease& operator=(const WorkspaceLease&) = delete;

  Workspace* get() const noexcept { return ws_.get(); }
  Workspace* operator->() const noexcept { return ws_.get(); }
  Workspace& operator*() const noexcept { return *ws_; }
  explicit operator bool() const noexcept { return ws_ != nullptr; }

  void reset() noexcept;

private:
  WorkspacePool* pool_ = nullptr;
  std::unique_ptr<Workspace> ws_;
};

class WorkspacePool {
public:
  explicit WorkspacePool(WorkspaceBudget budget, std::size_t max_capacity = 8);
  ~WorkspacePool() = default;

  WorkspacePool(const WorkspacePool&) = delete;
  WorkspacePool& operator=(const WorkspacePool&) = delete;

  WorkspaceLease acquire();

  const WorkspaceBudget& budget() const noexcept { return budget_; }
  std::size_t max_capacity() const noexcept { return max_capacity_; }
  std::size_t active_count() const noexcept;
  std::size_t idle_count() const noexcept;

private:
  friend class WorkspaceLease;
  void release(std::unique_ptr<Workspace> ws) noexcept;

  WorkspaceBudget budget_;
  std::size_t max_capacity_;
  mutable std::mutex mutex_;
  std::condition_variable cv_;
  std::vector<std::unique_ptr<Workspace>> idle_;
  std::size_t active_count_ = 0;
};

} // namespace neo_fft::runtime
