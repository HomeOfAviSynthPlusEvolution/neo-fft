#include "spectral/fft.hpp"
#include <pocketfft_hdronly.h>
#include <algorithm>

namespace neo_fft {
namespace {
template <class T> std::size_t validate(T* data, BatchLayout l, int h, int w) {
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
template <class T> void scan(const T* p, BatchLayout l, int h, int w) {
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
  for (std::size_t b = 0; b < l.active; ++b)
    for (int x : {0, w % 2 == 0 ? w / 2 : 0})
      for (int y = 0; y < h; ++y) {
        const auto a = p[b * l.distance + y * l.row_stride + x];
        const auto c = std::conj(p[b * l.distance + ((h - y) % h) * l.row_stride + x]);
        // FFT roundoff on boundary columns is allowed; arbitrary spectra are not.
        const double scale = std::max({1.0, std::abs(std::complex<double>(a)), std::abs(std::complex<double>(c))});
        require(std::abs(std::complex<double>(a) - std::complex<double>(c)) <= 64 * 1.1920928955078125e-7 * scale,
                "inverse FFT requires Hermitian-compatible boundary columns");
      }
}
} // namespace
RealFFT::RealFFT(int height, int width) : height_(dimension(height)), width_(dimension(width)) {
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
  const pocketfft::shape_t shape{std::size_t(height_), std::size_t(width_)}, axes{0, 1};
  const pocketfft::stride_t rs{std::ptrdiff_t(r.row_stride * sizeof(float)), sizeof(float)};
  const pocketfft::stride_t ss{std::ptrdiff_t(s.row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  for (std::size_t b = 0; b < r.active; ++b)
    pocketfft::r2c(shape, rs, ss, axes, true, in + b * r.distance, out + b * s.distance, 1.0f, 1);
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
  const pocketfft::shape_t shape{std::size_t(height_), std::size_t(width_)}, axes{0, 1};
  const pocketfft::stride_t rs{std::ptrdiff_t(r.row_stride * sizeof(float)), sizeof(float)};
  const pocketfft::stride_t ss{std::ptrdiff_t(s.row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  for (std::size_t b = 0; b < r.active; ++b)
    pocketfft::c2r(shape, ss, rs, axes, false, in + b * s.distance, out + b * r.distance,
                  1.0f / (float(width_) * float(height_)), 1);
  scan(out, r, height_, width_);
}
void RealFFT::forward(const float* in, std::complex<float>* out) const {
  forward(in, {std::size_t(width_), samples(), 1, 1}, out, {std::size_t(columns()), bins(), 1, 1});
}
void RealFFT::inverse(const std::complex<float>* in, float* out) const {
  inverse(in, {std::size_t(columns()), bins(), 1, 1}, out, {std::size_t(width_), samples(), 1, 1});
}
} // namespace neo_fft
