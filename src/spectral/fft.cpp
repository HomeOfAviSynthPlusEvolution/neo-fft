#include "spectral/fft.hpp"
#include "spectral/fft_backend.hpp"
#include <algorithm>

#if NEO_FFT_ENABLE_SIMD
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
#include <hwy/targets.h>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif

namespace neo_fft {
namespace detail {
const FftBackend& native_fft() noexcept {
#if NEO_FFT_ENABLE_SIMD
#if NEO_FFT_FFT_X86_TARGETS
  const auto targets = hwy::SupportedTargets();
  if (targets & (HWY_AVX3 | HWY_AVX3_SPR | HWY_AVX3_ZEN4 | HWY_AVX3_DL)) {
    if (avx512_fft().lanes > 1)
      return avx512_fft();
  }
  if (targets & HWY_AVX2) {
    if (avx2_fft().lanes > 1)
      return avx2_fft();
  }
  if (targets & (HWY_SSE2 | HWY_SSSE3 | HWY_SSE4)) {
    if (sse2_fft().lanes > 1)
      return sse2_fft();
  }
#else
  if (native_target_fft().lanes > 1)
    return native_target_fft();
#endif
#endif
  return scalar_fft();
}
} // namespace detail

const FftBackend& backend_by_profile(FftProfile profile) noexcept {
  switch (profile) {
    case FftProfile::scalar:
      return detail::scalar_fft();
#if NEO_FFT_FFT_X86_TARGETS
    case FftProfile::sse2:
      return detail::sse2_fft();
    case FftProfile::avx2:
      return detail::avx2_fft();
    case FftProfile::avx512:
      return detail::avx512_fft();
#endif
    case FftProfile::native:
    default:
      return detail::native_fft();
  }
}

int fft_lanes(FftProfile profile) noexcept {
  return backend_by_profile(profile).lanes;
}

const char* fft_profile_name(FftProfile profile) noexcept {
  switch (profile) {
    case FftProfile::scalar:
      return "pocketfft-scalar";
    case FftProfile::sse2:
      return "pocketfft-sse2";
    case FftProfile::avx2:
      return "pocketfft-avx2";
    case FftProfile::avx512:
      return "pocketfft-avx512";
    case FftProfile::native:
    default:
      return fft_lanes(profile) > 1 ? "pocketfft-native" : "pocketfft-scalar";
  }
}

const char* fft_backend_name(FftProfile profile) noexcept {
  return backend_by_profile(profile).name;
}

bool fft_profile_supported(FftProfile profile) noexcept {
  switch (profile) {
    case FftProfile::scalar:
      return true;
#if NEO_FFT_ENABLE_SIMD
#if NEO_FFT_FFT_X86_TARGETS
    case FftProfile::sse2:
      return (hwy::SupportedTargets() & (HWY_SSE2 | HWY_SSSE3 | HWY_SSE4)) != 0;
    case FftProfile::avx2:
      return (hwy::SupportedTargets() & HWY_AVX2) != 0;
    case FftProfile::avx512:
      return (hwy::SupportedTargets() & (HWY_AVX3 | HWY_AVX3_SPR | HWY_AVX3_ZEN4 | HWY_AVX3_DL)) != 0;
#endif
    case FftProfile::native:
      return true;
#endif
    default:
      return profile == FftProfile::scalar || profile == FftProfile::native;
  }
}

RealFFT::RealFFT(int height, int width, FftProfile profile)
    : height_(dimension(height)), width_(dimension(width)), profile_(profile), backend_(backend_by_profile(profile)) {
  require(fft_profile_supported(profile), "FFT profile is not supported by current CPU");
  plane_extent<float>(width_, height_, static_cast<std::ptrdiff_t>(mul_size(width_, sizeof(float))));
  plane_extent<std::complex<float>>(columns(), height_,
                                    static_cast<std::ptrdiff_t>(mul_size(columns(), sizeof(std::complex<float>))));
}
void RealFFT::forward(const float* in, BatchLayout r, std::complex<float>* out, BatchLayout s) const {
  require(r.active == s.active, "FFT batch counts differ");
  if (!r.active)
    return;
  backend_.batch_r2c(r.active, height_, width_, in, r.distance, r.row_stride, out, s.distance, s.row_stride);
}
void RealFFT::inverse(const std::complex<float>* in, BatchLayout s, float* out, BatchLayout r) const {
  require(r.active == s.active, "FFT batch counts differ");
  if (!r.active)
    return;
  const float scale = 1.0f / (float(width_) * float(height_));
  backend_.batch_c2r(r.active, height_, width_, in, s.distance, s.row_stride, out, r.distance, r.row_stride, scale);
}
void RealFFT::forward(const float* in, std::complex<float>* out) const {
  forward(in, {std::size_t(width_), samples(), 1, 1}, out, {std::size_t(columns()), bins(), 1, 1});
}
void RealFFT::inverse(const std::complex<float>* in, float* out) const {
  inverse(in, {std::size_t(columns()), bins(), 1, 1}, out, {std::size_t(width_), samples(), 1, 1});
}
} // namespace neo_fft
