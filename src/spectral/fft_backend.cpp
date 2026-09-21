#include "spectral/fft_backend.hpp"

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
#include <pocketfft_hdronly.h>

namespace neo_fft::detail {
namespace {
namespace pf = POCKETFFT_NAMESPACE;

void r2c(int height, int width, const float* in, std::size_t in_row_stride, std::complex<float>* out,
         std::size_t out_row_stride) {
  const pf::shape_t shape{std::size_t(height), std::size_t(width)}, axes{0, 1};
  const pf::stride_t rs{std::ptrdiff_t(in_row_stride * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void c2r(int height, int width, const std::complex<float>* in, std::size_t in_row_stride, float* out,
         std::size_t out_row_stride, float fct) {
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
  if (batch == 1) {
    r2c(height, width, in, in_row_stride, out, out_row_stride);
    return;
  }
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
  if (batch == 1) {
    c2r(height, width, in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  const pf::shape_t shape{batch, std::size_t(height), std::size_t(width)}, axes{1, 2};
  const pf::stride_t rs{std::ptrdiff_t(out_dist * sizeof(float)), std::ptrdiff_t(out_row_stride * sizeof(float)),
                        sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(in_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}
} // namespace

const FftBackend& BACKEND_FN() noexcept {
  static const FftBackend backend{int(pf::detail::VLEN<float>::val), BACKEND_NAME, r2c, c2r, batch_r2c, batch_c2r};
  return backend;
}

} // namespace neo_fft::detail
