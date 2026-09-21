#pragma once

#include <complex>
#include <cstddef>

namespace neo_fft::codelet {

using cfloat = std::complex<float>;

#if defined(NEO_FFT_HAS_AVX2_CODELET)

void batch_fft8x8_r2c(std::size_t batch, const float* in, std::size_t in_dist,
                      std::size_t in_stride, cfloat* out, std::size_t out_dist,
                      std::size_t out_stride) noexcept;

void batch_fft8x8_c2r(std::size_t batch, const cfloat* in, std::size_t in_dist,
                      std::size_t in_stride, float* out, std::size_t out_dist,
                      std::size_t out_stride, float fct) noexcept;

inline void fft8x8_r2c(const float* in, std::size_t in_stride, cfloat* out,
                       std::size_t out_stride) noexcept {
  batch_fft8x8_r2c(1, in, 64, in_stride, out, 40, out_stride);
}

inline void fft8x8_c2r(const cfloat* in, std::size_t in_stride, float* out,
                       std::size_t out_stride, float fct) noexcept {
  batch_fft8x8_c2r(1, in, 40, in_stride, out, 64, out_stride, fct);
}

#endif

} // namespace neo_fft::codelet
