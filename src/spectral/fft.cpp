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

namespace {
template <class T>
std::size_t validate(T* data, BatchLayout l, int h, int w) {
  require(l.active <= l.capacity, "FFT active count exceeds capacity");
  if (!l.active)
    return 0;
  require(l.row_stride <= INT32_MAX, "FFT row stride exceeds int32");
  const auto stride = mul_size(l.row_stride, sizeof(T));
  require(stride <= PTRDIFF_MAX, "FFT row stride overflow");
  const auto slice = plane_extent<T>(w, h, static_cast<std::ptrdiff_t>(stride));
  const auto distance = mul_size(l.distance, sizeof(T));
  require(distance >= slice, "FFT batch distance shorter than slice");
  const auto total = add_size(mul_size(l.capacity - 1, distance), slice);
  require(total <= PTRDIFF_MAX, "FFT batch extent overflow");
  checked_plane(data, w, h, static_cast<std::ptrdiff_t>(stride), slice);
  require(total <= UINTPTR_MAX - reinterpret_cast<std::uintptr_t>(data), "FFT address extent overflow");
  return add_size(mul_size(l.active - 1, distance), slice);
}
template <class T>
void scan(const T* p, BatchLayout l, int h, int w) {
  for (std::size_t b = 0; b < l.active; ++b)
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x) {
        const auto v = p[b * l.distance + y * l.row_stride + x];
        if constexpr (std::is_same_v<T, float>) {
          finite(v);
        } else {
          finite(v.real());
          finite(v.imag());
        }
      }
}
void hermitian(const std::complex<float>* p, BatchLayout l, int h, int w) {
  for (std::size_t b = 0; b < l.active; ++b) {
    double scale = 1;
    for (int y = 0; y < h; ++y)
      for (int x = 0; x <= w / 2; ++x)
        scale = std::max(scale, std::abs(std::complex<double>(p[b * l.distance + y * l.row_stride + x])));
    for (int x : {0, w % 2 == 0 ? w / 2 : 0})
      for (int y = 0; y < h; ++y) {
        const auto a = p[b * l.distance + y * l.row_stride + x];
        const auto c = std::conj(p[b * l.distance + ((h - y) % h) * l.row_stride + x]);
        // FFT roundoff on boundary columns is allowed; arbitrary spectra are not.
        require(std::abs(std::complex<double>(a) - std::complex<double>(c)) <= 64 * 1.1920928955078125e-7 * scale,
                "inverse FFT requires Hermitian-compatible boundary columns");
      }
  }
}
} // namespace

RealFFT::RealFFT(int height, int width, FftProfile profile)
    : height_(dimension(height)), width_(dimension(width)), profile_(profile), backend_(backend_by_profile(profile)) {
  plane_extent<float>(width_, height_, static_cast<std::ptrdiff_t>(mul_size(width_, sizeof(float))));
  plane_extent<std::complex<float>>(columns(), height_,
                                    static_cast<std::ptrdiff_t>(mul_size(columns(), sizeof(std::complex<float>))));
}
void RealFFT::forward(const float* in, BatchLayout r, std::complex<float>* out, BatchLayout s) const {
  require(r.active == s.active, "FFT batch counts differ");
  const auto rn = validate(in, r, height_, width_), sn = validate(out, s, height_, columns());
  if (!r.active)
    return;
  disjoint(in, rn, out, sn);
  scan(in, r, height_, width_);
  for (std::size_t b = 0; b < r.active; ++b)
    backend_.r2c(height_, width_, in + b * r.distance, r.row_stride, out + b * s.distance, s.row_stride);
  scan(out, s, height_, columns());
}
void RealFFT::inverse(const std::complex<float>* in, BatchLayout s, float* out, BatchLayout r) const {
  require(r.active == s.active, "FFT batch counts differ");
  const auto rn = validate(out, r, height_, width_), sn = validate(in, s, height_, columns());
  if (!r.active)
    return;
  disjoint(in, sn, out, rn);
  scan(in, s, height_, columns());
  hermitian(in, s, height_, width_);
  const float scale = 1.0f / (float(width_) * float(height_));
  for (std::size_t b = 0; b < r.active; ++b)
    backend_.c2r(height_, width_, in + b * s.distance, s.row_stride, out + b * r.distance, r.row_stride, scale);
  scan(out, r, height_, width_);
}
void RealFFT::forward(const float* in, std::complex<float>* out) const {
  forward(in, {std::size_t(width_), samples(), 1, 1}, out, {std::size_t(columns()), bins(), 1, 1});
}
void RealFFT::inverse(const std::complex<float>* in, float* out) const {
  inverse(in, {std::size_t(columns()), bins(), 1, 1}, out, {std::size_t(width_), samples(), 1, 1});
}
} // namespace neo_fft
