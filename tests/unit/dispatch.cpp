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
    Guarded huge(33);
    for (int i = 0; i < 33; ++i)
      huge.data[i] = {1e30f, 0};
    rejects([&] { optimized(huge.data, nullptr, 33, 0, {}); });

    Guarded overflow_test(33);
    for (int i = 0; i < 33; ++i)
      overflow_test.data[i] = {2.0f, 0};
    SpectralParams overflow_params;
    overflow_params.type = 2;
    overflow_params.a = 3e38f;
    rejects([&] { optimized(overflow_test.data, nullptr, 33, 0, overflow_params); });
    std::cout << "spectral target: " << spectral_target(0) << "; scalar comparison and protected-page tails passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
