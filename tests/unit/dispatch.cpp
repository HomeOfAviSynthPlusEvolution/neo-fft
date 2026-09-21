#include "kernels/spectral.hpp"
#include "../test.hpp"
#include <algorithm>
#include <cstring>
#include <memory>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
using namespace neo_fft;
struct Guarded {
  std::complex<float>* data = nullptr;
  void* allocation = nullptr;
  std::vector<std::complex<float>> fallback;
  explicit Guarded(std::size_t count) {
#if defined(_WIN32)
    SYSTEM_INFO info{};
    GetSystemInfo(&info);
    CHECK(count * sizeof(std::complex<float>) <= info.dwPageSize);
    allocation = VirtualAlloc(nullptr, 2 * info.dwPageSize, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    CHECK(allocation);
    auto* boundary = static_cast<unsigned char*>(allocation) + info.dwPageSize;
    DWORD previous = 0;
    CHECK(VirtualProtect(boundary, info.dwPageSize, PAGE_NOACCESS, &previous));
    data = reinterpret_cast<std::complex<float>*>(boundary - count * sizeof(std::complex<float>));
    for (std::size_t i = 0; i < count; ++i)
      new (data + i) std::complex<float>();
#else
    fallback.resize(count);
    data = fallback.data();
#endif
  }
  ~Guarded() {
#if defined(_WIN32)
    if (allocation)
      VirtualFree(allocation, 0, MEM_RELEASE);
#endif
  }
  Guarded(const Guarded&) = delete;
  Guarded& operator=(const Guarded&) = delete;
};
int main() {
  try {
    const auto optimized = select_spectral(0);
    optimized(nullptr, nullptr, 0, 0, {});
    for (int type = -1; type <= 4; ++type)
      for (float exponent : {.5f, 1.0f, 2.0f, .50005f})
        for (std::size_t count : {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 257})
          for (bool mean : {false, true}) {
            Guarded samples(count), grid(count);
            auto expected = buffer<std::complex<float>>(count), original_grid = expected;
            for (std::size_t i = 0; i < count; ++i) {
              samples.data[i] = expected[i] = {float(int(i % 19) - 9), float(int(i % 23) - 11)};
              grid.data[i] = original_grid[i] = {float(i % 5) * .25f, float(i % 3) * .5f};
            }
            SpectralParams p{type, 9, .5f, 25, 100, exponent, .25f};
            spectral_scalar(expected.data(), mean ? grid.data : nullptr, count, .5f, p);
            optimized(samples.data, mean ? grid.data : nullptr, count, .5f, p);
            CHECK(std::memcmp(grid.data, original_grid.data(), count * sizeof(std::complex<float>)) == 0);
            for (std::size_t i = 0; i < count; ++i) {
              check_near(samples.data[i].real(), expected[i].real(), 2e-6);
              check_near(samples.data[i].imag(), expected[i].imag(), 2e-6);
            }
          }
    // 1. Pure SIMD (exact vector length, count = 32): verify SIMD detection without scalar tail
    {
      SpectralParams overflow_params;
      overflow_params.type = 2;
      overflow_params.a = 3e38f;

      // grid == nullptr: pure SIMD, index 0 overflows (2 * 3e38 = Inf)
      Guarded pure_simd(32);
      pure_simd.data[0] = {2.0f, 0};
      rejects([&] { optimized(pure_simd.data, nullptr, 32, 0, overflow_params); });

      // grid != nullptr: pure SIMD, index 0 overflows (re * gain + mr = Inf)
      Guarded pure_simd_grid(32), grid(32);
      pure_simd_grid.data[0] = {2.0f, 0};
      rejects([&] { optimized(pure_simd_grid.data, grid.data, 32, 0.5f, overflow_params); });
    }

    // 2. SIMD overflow with scalar tail present (count = 33): tail elements are finite & non-overflowing
    {
      SpectralParams overflow_params;
      overflow_params.type = 2;
      overflow_params.a = 3e38f;

      // grid == nullptr: only index 0 overflows, index 32 (tail) has 0.0f (gain * 0 = 0, finite)
      Guarded overflow_test(33);
      overflow_test.data[0] = {2.0f, 0};
      rejects([&] { optimized(overflow_test.data, nullptr, 33, 0, overflow_params); });

      // grid != nullptr: only index 0 overflows in SIMD, index 32 (tail) is finite
      Guarded overflow_test_grid(33), grid(33);
      overflow_test_grid.data[0] = {2.0f, 0};
      rejects([&] { optimized(overflow_test_grid.data, grid.data, 33, 0.5f, overflow_params); });
    }

    // 3. Huge power overflow: only index 0 has huge amplitude, scalar tail is normal and finite
    {
      Guarded huge_no_grid(33);
      huge_no_grid.data[0] = {1e30f, 0};
      for (std::size_t i = 1; i < 33; ++i)
        huge_no_grid.data[i] = {1.0f, 0};
      rejects([&] { optimized(huge_no_grid.data, nullptr, 33, 0, {}); });

      Guarded huge_grid(33), grid(33);
      huge_grid.data[0] = {1e30f, 0};
      for (std::size_t i = 1; i < 33; ++i)
        huge_grid.data[i] = {1.0f, 0};
      rejects([&] { optimized(huge_grid.data, grid.data, 33, 0.5f, {}); });
    }
    std::cout << "spectral target: " << spectral_target(0) << "; scalar comparison and protected-page tails passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
