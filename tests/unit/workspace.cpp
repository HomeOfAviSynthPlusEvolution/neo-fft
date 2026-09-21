#include "runtime/workspace.hpp"
#include "runtime/workspace_pool.hpp"
#include "../test.hpp"
#include <atomic>
#include <chrono>
#include <future>
#include <thread>
#include <vector>

using namespace neo_fft;
using namespace neo_fft::runtime;

int main() {
  try {
    // 1. Geometry and Budget verification
    Axis ax = fft3d_axis(1280, 32, 16);
    Axis ay = fft3d_axis(720, 32, 16);
    Geometry geom(ax, ay);
    RealFFT fft(32, 32);

    auto budget_fft3d = make_workspace_budget(geom, fft, true, 1);
    CHECK(budget_fft3d.accum_width == ax.cover);
    CHECK(budget_fft3d.accum_height == ay.cover);
    CHECK(budget_fft3d.accum_stride_bytes % kSimdAlignment == 0);
    CHECK(budget_fft3d.accum_offset % kSimdAlignment == 0);
    CHECK(budget_fft3d.row_width == ax.cover);
    CHECK(budget_fft3d.row_height == ay.block);
    CHECK(budget_fft3d.row_stride_bytes % kSimdAlignment == 0);
    CHECK(budget_fft3d.row_offset % kSimdAlignment == 0);
    CHECK(budget_fft3d.block_offset % kSimdAlignment == 0);
    CHECK(budget_fft3d.inverse_offset % kSimdAlignment == 0);
    CHECK(budget_fft3d.spectrum_offset % kSimdAlignment == 0);
    CHECK(budget_fft3d.total_bytes > 0);

    auto budget_dft = make_workspace_budget(geom, fft, false, 1);
    CHECK(budget_dft.row_height == 0);
    CHECK(budget_dft.row_stride_bytes == 0);

    // 2. Workspace SIMD 64-byte alignment and span2d views
    Workspace ws(budget_fft3d);
    CHECK(reinterpret_cast<std::uintptr_t>(ws.accum().data()) % kSimdAlignment == 0);
    CHECK(reinterpret_cast<std::uintptr_t>(ws.accum().row_ptr(1)) % kSimdAlignment == 0);
    CHECK(reinterpret_cast<std::uintptr_t>(ws.row().data()) % kSimdAlignment == 0);
    CHECK(reinterpret_cast<std::uintptr_t>(ws.block().data()) % kSimdAlignment == 0);
    CHECK(reinterpret_cast<std::uintptr_t>(ws.inverse().data()) % kSimdAlignment == 0);
    CHECK(reinterpret_cast<std::uintptr_t>(ws.spectrum().data()) % kSimdAlignment == 0);

    // span2d 2D planes and 1D spans
    auto acc = ws.accum();
    CHECK(acc.width() == ax.cover);
    CHECK(acc.height() == ay.cover);
    acc(0, 0) = 42.0f;
    acc(5, 10) = 99.0f;
    CHECK(acc(0, 0) == 42.0f);
    CHECK(acc(5, 10) == 99.0f);

    // Subplane view test
    auto sub = acc.subplane(10, 5, 8, 8);
    CHECK(sub.width() == 8 && sub.height() == 8);
    CHECK(sub(0, 0) == 99.0f);
    sub(1, 1) = 123.0f;
    CHECK(acc(6, 11) == 123.0f);

    // Restrict views
    auto acc_restr = ws.accum_restrict();
    CHECK(acc_restr(0, 0) == 42.0f);

    // Reset verification (zero logical extent)
    ws.reset();
    CHECK(acc(0, 0) == 0.0f);
    CHECK(acc(5, 10) == 0.0f);
    CHECK(acc(6, 11) == 0.0f);

    // 3. WorkspacePool and RAII WorkspaceLease
    WorkspacePool pool(budget_fft3d, 4);
    CHECK(pool.max_capacity() == 4);
    CHECK(pool.active_count() == 0);
    CHECK(pool.idle_count() == 0);

    {
      auto lease = pool.acquire();
      CHECK(bool(lease));
      CHECK(pool.active_count() == 1);
      CHECK(pool.idle_count() == 0);
      lease->accum()(0, 0) = 3.14f;

      // Move construction
      WorkspaceLease lease2 = std::move(lease);
      CHECK(!bool(lease));
      CHECK(bool(lease2));
      CHECK(lease2->accum()(0, 0) == 3.14f);
      CHECK(pool.active_count() == 1);
    }
    // Automatically released to pool
    CHECK(pool.active_count() == 0);
    CHECK(pool.idle_count() == 1);

    // Reuse cached workspace
    {
      auto lease = pool.acquire();
      CHECK(pool.active_count() == 1);
      CHECK(pool.idle_count() == 0);
      // Value before reset
      CHECK(lease->accum()(0, 0) == 3.14f);
      lease->reset();
      CHECK(lease->accum()(0, 0) == 0.0f);
    }
    CHECK(pool.active_count() == 0);
    CHECK(pool.idle_count() == 1);

    // 4. Concurrent stress test
    constexpr int kThreads = 8;
    constexpr int kIters = 25;
    std::vector<std::future<void>> futures;
    for (int t = 0; t < kThreads; ++t) {
      futures.push_back(std::async(std::launch::async, [&pool, t] {
        for (int i = 0; i < kIters; ++i) {
          auto lease = pool.acquire();
          lease->reset();
          auto a = lease->accum();
          const float val = float(t * 1000 + i + 1);
          a(0, 0) = val;
          a(10, 10) = val * 2.0f;
          std::this_thread::yield();
          CHECK(a(0, 0) == val);
          CHECK(a(10, 10) == val * 2.0f);
        }
      }));
    }
    for (auto& f : futures) {
      f.get();
    }
    CHECK(pool.active_count() == 0);
    CHECK(pool.idle_count() <= pool.max_capacity());

    // 5. Max capacity blocking test
    WorkspacePool small_pool(budget_dft, 2);
    auto l1 = small_pool.acquire();
    auto l2 = small_pool.acquire();
    CHECK(small_pool.active_count() == 2);

    std::atomic<bool> acquired_l3{false};
    auto future_l3 = std::async(std::launch::async, [&small_pool, &acquired_l3] {
      auto l3 = small_pool.acquire();
      acquired_l3.store(true);
      return bool(l3);
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    CHECK(acquired_l3.load() == false);

    // Release l1 to allow l3 to acquire
    l1.reset();
    CHECK(future_l3.get() == true);
    CHECK(acquired_l3.load() == true);

    std::cout << "workspace unit tests passed successfully.\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "FAILED: " << e.what() << "\n";
    return 1;
  }
}
