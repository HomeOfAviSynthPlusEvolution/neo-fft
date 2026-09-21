#include "spectral/fft_backend.hpp"
#include "spectral/codelets/fft16x16.hpp"
#include "spectral/codelets/fft8x8.hpp"
#include "spectral/codelets/fft32x32.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
class TlsScratchPool {
public:
  static constexpr std::size_t CAPACITY = 1024 * 1024; // 1 MB

  TlsScratchPool() = default;
  ~TlsScratchPool() {
    if (buffer_) {
#if defined(_WIN32)
      _aligned_free(buffer_);
#else
      std::free(buffer_);
#endif
      buffer_ = nullptr;
    }
  }

  TlsScratchPool(const TlsScratchPool&) = delete;
  TlsScratchPool& operator=(const TlsScratchPool&) = delete;

  void* allocate(std::size_t size) {
    if (active_) {
      if (!buffer_) {
#if defined(_WIN32)
        buffer_ = static_cast<char*>(_aligned_malloc(CAPACITY, 64));
#else
        void* ptr = nullptr;
        if (posix_memalign(&ptr, 64, CAPACITY) != 0) {
          ptr = std::malloc(CAPACITY);
        }
        buffer_ = static_cast<char*>(ptr);
#endif
      }
      if (buffer_) {
        const std::size_t aligned_size = (size + 63) & ~std::size_t(63);
        if (offset_ + aligned_size <= CAPACITY) {
          char* ptr = buffer_ + offset_;
          offset_ += aligned_size;
          ++active_allocs_;
          last_ptr_ = ptr;
          last_size_ = aligned_size;
          return ptr;
        }
      }
    }
    return std::malloc(size);
  }

  void deallocate(void* ptr) noexcept {
    if (!ptr)
      return;
    if (buffer_ && ptr >= buffer_ && ptr < buffer_ + CAPACITY) {
      if (active_allocs_ > 0) {
        --active_allocs_;
      }
      if (active_allocs_ == 0) {
        offset_ = 0;
        last_ptr_ = nullptr;
        last_size_ = 0;
      } else if (ptr == last_ptr_) {
        offset_ -= last_size_;
        last_ptr_ = nullptr;
        last_size_ = 0;
      }
      return;
    }
    std::free(ptr);
  }

  bool is_active() const noexcept {
    return active_;
  }

  void set_active(bool active) noexcept {
    active_ = active;
  }

  std::size_t active_allocs() const noexcept {
    return active_allocs_;
  }

  void reset() noexcept {
    if (active_allocs_ == 0) {
      offset_ = 0;
      last_ptr_ = nullptr;
      last_size_ = 0;
    }
  }

private:
  char* buffer_ = nullptr;
  std::size_t offset_ = 0;
  std::size_t active_allocs_ = 0;
  char* last_ptr_ = nullptr;
  std::size_t last_size_ = 0;
  bool active_ = false;
};

thread_local TlsScratchPool g_scratch_pool;

inline void* pocketfft_scratch_malloc(std::size_t size) {
  return g_scratch_pool.allocate(size);
}

inline void pocketfft_scratch_free(void* ptr) noexcept {
  g_scratch_pool.deallocate(ptr);
}

class ScratchScope {
public:
  ScratchScope() {
    g_scratch_pool.set_active(true);
  }
  ~ScratchScope() {
    g_scratch_pool.set_active(false);
    g_scratch_pool.reset();
  }
  ScratchScope(const ScratchScope&) = delete;
  ScratchScope& operator=(const ScratchScope&) = delete;
};

class SuspendScratchScope {
public:
  SuspendScratchScope() noexcept : prev_(g_scratch_pool.is_active()) {
    g_scratch_pool.set_active(false);
  }
  ~SuspendScratchScope() noexcept {
    g_scratch_pool.set_active(prev_);
  }
  SuspendScratchScope(const SuspendScratchScope&) = delete;
  SuspendScratchScope& operator=(const SuspendScratchScope&) = delete;

private:
  bool prev_;
};
} // namespace

// Built multiple times with separate dependency namespaces. No shared
// cache or worker pool exists between targets.
#if defined(NEO_FFT_FFT_SCALAR)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_scalar_c90e55b3
#ifndef POCKETFFT_NO_VECTORS
#define POCKETFFT_NO_VECTORS
#endif
#define BACKEND_FN scalar_fft
#define BACKEND_NAME "pocketfft-scalar"
#elif defined(NEO_FFT_FFT_SSE2)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_sse2_c90e55b3
#define BACKEND_FN sse2_fft
#define BACKEND_NAME "pocketfft-sse2"
#elif defined(NEO_FFT_FFT_AVX2)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_avx2_c90e55b3
#define BACKEND_FN avx2_fft
#define BACKEND_NAME "pocketfft-avx2"
#elif defined(NEO_FFT_FFT_AVX512)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_avx512_c90e55b3
#define BACKEND_FN avx512_fft
#define BACKEND_NAME "pocketfft-avx512"
#else
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_native_c90e55b3
#define BACKEND_FN native_target_fft
#define BACKEND_NAME "pocketfft-native"
#endif
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_CACHE_SIZE 16

#define malloc(sz) pocketfft_scratch_malloc(sz)
#define free(ptr) pocketfft_scratch_free(ptr)
#include <pocketfft_hdronly.h>
#undef malloc
#undef free

namespace POCKETFFT_NAMESPACE::detail {

template <typename T>
inline std::shared_ptr<T> get_plan_persistent(size_t length) {
  static std::mutex mut;
  static std::unordered_map<size_t, std::shared_ptr<T>> cache;
  {
    std::lock_guard<std::mutex> lock(mut);
    auto it = cache.find(length);
    if (it != cache.end()) {
      return it->second;
    }
  }
  SuspendScratchScope suspend;
  auto plan = std::make_shared<T>(length);
  {
    std::lock_guard<std::mutex> lock(mut);
    cache.emplace(length, plan);
  }
  return plan;
}

template <>
inline std::shared_ptr<pocketfft_r<float>> get_plan<pocketfft_r<float>>(size_t length) {
  return get_plan_persistent<pocketfft_r<float>>(length);
}

template <>
inline std::shared_ptr<pocketfft_c<float>> get_plan<pocketfft_c<float>>(size_t length) {
  return get_plan_persistent<pocketfft_c<float>>(length);
}

template <>
inline std::shared_ptr<pocketfft_r<double>> get_plan<pocketfft_r<double>>(size_t length) {
  return get_plan_persistent<pocketfft_r<double>>(length);
}

template <>
inline std::shared_ptr<pocketfft_c<double>> get_plan<pocketfft_c<double>>(size_t length) {
  return get_plan_persistent<pocketfft_c<double>>(length);
}

} // namespace POCKETFFT_NAMESPACE::detail

namespace neo_fft::detail {
namespace {
namespace pf = POCKETFFT_NAMESPACE;

inline void ensure_plans(int height, int width) {
  pf::detail::get_plan<pf::detail::pocketfft_r<float>>(std::size_t(width));
  pf::detail::get_plan<pf::detail::pocketfft_r<float>>(std::size_t(height));
  pf::detail::get_plan<pf::detail::pocketfft_c<float>>(std::size_t(width));
  pf::detail::get_plan<pf::detail::pocketfft_c<float>>(std::size_t(height));
}

void r2c(int height, int width, const float* in, std::size_t in_row_stride, std::complex<float>* out,
         std::size_t out_row_stride) {
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::fft32x32_r2c(in, in_row_stride, out, out_row_stride);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::fft16x16_r2c(in, in_row_stride, out, out_row_stride);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::fft8x8_r2c(in, in_row_stride, out, out_row_stride);
    return;
  }
#endif
  ensure_plans(height, width);
  ScratchScope scope;
  const pf::shape_t shape{std::size_t(height), std::size_t(width)}, axes{0, 1};
  const pf::stride_t rs{std::ptrdiff_t(in_row_stride * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void c2r(int height, int width, const std::complex<float>* in, std::size_t in_row_stride, float* out,
         std::size_t out_row_stride, float fct) {
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::fft32x32_c2r(in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::fft16x16_c2r(in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::fft8x8_c2r(in, in_row_stride, out, out_row_stride, fct);
    return;
  }
#endif
  ensure_plans(height, width);
  ScratchScope scope;
  const pf::shape_t shape{std::size_t(height), std::size_t(width)}, axes{0, 1};
  const pf::stride_t rs{std::ptrdiff_t(out_row_stride * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}

void batch_r2c(std::size_t batch, int height, int width, const float* in, std::size_t in_dist,
               std::size_t in_row_stride, std::complex<float>* out, std::size_t out_dist,
               std::size_t out_row_stride) {
  if (batch == 0)
    return;
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::batch_fft32x32_r2c(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::batch_fft16x16_r2c(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::batch_fft8x8_r2c(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride);
    return;
  }
#endif
  if (batch == 1) {
    r2c(height, width, in, in_row_stride, out, out_row_stride);
    return;
  }
  ensure_plans(height, width);
  ScratchScope scope;
  const pf::shape_t shape{batch, std::size_t(height), std::size_t(width)}, axes{1, 2};
  const pf::stride_t rs{std::ptrdiff_t(in_dist * sizeof(float)), std::ptrdiff_t(in_row_stride * sizeof(float)),
                        sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(out_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void batch_c2r(std::size_t batch, int height, int width, const std::complex<float>* in, std::size_t in_dist,
               std::size_t in_row_stride, float* out, std::size_t out_dist, std::size_t out_row_stride, float fct) {
  if (batch == 0)
    return;
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::batch_fft32x32_c2r(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride, fct);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::batch_fft16x16_c2r(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride, fct);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::batch_fft8x8_c2r(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride, fct);
    return;
  }
#endif
  if (batch == 1) {
    c2r(height, width, in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  ensure_plans(height, width);
  ScratchScope scope;
  const pf::shape_t shape{batch, std::size_t(height), std::size_t(width)}, axes{1, 2};
  const pf::stride_t rs{std::ptrdiff_t(out_dist * sizeof(float)), std::ptrdiff_t(out_row_stride * sizeof(float)),
                        sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(in_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}

inline void ensure_plans_3d(int depth, int height, int width) {
  ensure_plans(height, width);
  pf::detail::get_plan<pf::detail::pocketfft_r<float>>(std::size_t(depth));
  pf::detail::get_plan<pf::detail::pocketfft_c<float>>(std::size_t(depth));
}

void r2c_3d(int depth, int height, int width, const float* in, std::complex<float>* out) {
  ensure_plans_3d(depth, height, width);
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{0, 1, 2};
  const pf::stride_t rs{std::ptrdiff_t(height * width * sizeof(float)), std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(height * k * sizeof(std::complex<float>)), std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void c2r_3d(int depth, int height, int width, const std::complex<float>* in, float* out, float fct) {
  ensure_plans_3d(depth, height, width);
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{0, 1, 2};
  const pf::stride_t rs{std::ptrdiff_t(height * width * sizeof(float)), std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(height * k * sizeof(std::complex<float>)), std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}

void batch_r2c_3d(std::size_t batch, int depth, int height, int width, const float* in, std::size_t in_dist,
                  std::complex<float>* out, std::size_t out_dist) {
  if (batch == 0)
    return;
  if (batch == 1) {
    r2c_3d(depth, height, width, in, out);
    return;
  }
  ensure_plans_3d(depth, height, width);
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{batch, std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{1, 2, 3};
  const pf::stride_t rs{std::ptrdiff_t(in_dist * sizeof(float)), std::ptrdiff_t(height * width * sizeof(float)),
                        std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(height * k * sizeof(std::complex<float>)),
                        std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void batch_c2r_3d(std::size_t batch, int depth, int height, int width, const std::complex<float>* in, std::size_t in_dist,
                  float* out, std::size_t out_dist, float fct) {
  if (batch == 0)
    return;
  if (batch == 1) {
    c2r_3d(depth, height, width, in, out, fct);
    return;
  }
  ensure_plans_3d(depth, height, width);
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{batch, std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{1, 2, 3};
  const pf::stride_t rs{std::ptrdiff_t(out_dist * sizeof(float)), std::ptrdiff_t(height * width * sizeof(float)),
                        std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(height * k * sizeof(std::complex<float>)),
                        std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}
} // namespace

const FftBackend& BACKEND_FN() noexcept {
  static const FftBackend backend{int(pf::detail::VLEN<float>::val), BACKEND_NAME, r2c, c2r, batch_r2c, batch_c2r,
                                  r2c_3d, c2r_3d, batch_r2c_3d, batch_c2r_3d};
  return backend;
}

} // namespace neo_fft::detail
