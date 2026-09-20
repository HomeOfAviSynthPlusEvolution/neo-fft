#include "algorithms/plan.hpp"
#include <type_traits>

namespace neo_fft {
namespace {
void nonnegative(float v, const char* name) {
  require(std::isfinite(v) && v >= 0, std::string("invalid ") + name);
}
void valid_format(SampleFormat f) {
  require(f.floating ? f.bits == 32 : (f.bits == 8 || f.bits == 10 || f.bits == 12 || f.bits == 14 || f.bits == 16),
          "unsupported sample format");
}
Geometry geometry3d(int w, int h, const FFT3DConfig& c) {
  validate(c);
  return {fft3d_axis(w, c.bw, c.ow), fft3d_axis(h, c.bh, c.oh)};
}
Geometry geometry_dft(int w, int h, const DFTConfig& c) {
  validate(c);
  return {dft_axis(w, c.block, c.mode, c.overlap), dft_axis(h, c.block, c.mode, c.overlap)};
}
} // namespace
void validate(const FFT3DConfig& c) {
  require(c.bw >= 2 && c.bh >= 2, "FFT3D bw/bh must be >=2");
  require(c.ow <= c.bw / 2 && c.oh <= c.bh / 2, "FFT3D ow/oh exceeds half block");
  require(c.wintype >= 0 && c.wintype <= 2, "FFT3D wintype outside 0..2");
  nonnegative(c.sigma, "FFT3D sigma");
  nonnegative(c.degrid, "FFT3D degrid");
  require(std::isfinite(c.beta) && c.beta >= 1, "FFT3D beta must be >=1");
  select_spectral(c.opt);
}
void validate(const DFTConfig& c) {
  require(c.block > 0 && (c.mode == 0 || c.mode == 1), "DFTTest invalid sbsize/smode");
  require(c.mode != 0 || c.block % 2 == 1, "DFTTest center mode requires odd sbsize");
  require(c.mode == 0 || (c.overlap >= 0 && c.overlap < c.block), "DFTTest invalid sosize");
  if (c.mode == 1 && c.overlap > c.block / 2)
    require(c.block % (c.block - c.overlap) == 0, "DFTTest heavy overlap requires divisible step");
  require(c.swin >= 0 && c.swin <= 11 && c.twin >= 0 && c.twin <= 11, "DFTTest invalid swin/twin");
  require(c.ftype >= 0 && c.ftype <= 4, "DFTTest ftype outside 0..4");
  nonnegative(c.sbeta, "DFTTest sbeta");
  nonnegative(c.tbeta, "DFTTest tbeta");
  nonnegative(c.sigma, "DFTTest sigma");
  nonnegative(c.sigma2, "DFTTest sigma2");
  nonnegative(c.pmin, "DFTTest pmin");
  nonnegative(c.pmax, "DFTTest pmax");
  require(c.pmin <= c.pmax, "DFTTest pmin exceeds pmax");
  require(std::isfinite(c.f0beta) && c.f0beta > 0, "DFTTest f0beta must be positive");
  select_spectral(c.opt);
}
Plan::Plan(int w, int h, SampleFormat f, const FFT3DConfig& c)
    : geometry(geometry3d(w, h, c)), format(f), algorithm(Algorithm::FFT3D), fft(c.bh, c.bw),
      wx_(fft3d_window(c.bw, geometry.x.overlap, c.wintype)), wy_(fft3d_window(c.bh, geometry.y.overlap, c.wintype)),
      kernel_(select_spectral(c.opt)), mean_scale_(c.degrid) {
  valid_format(f);
  const float factor = f.floating ? 1.0f / 255 : float(1 << (f.bits - 8));
  const float sigma = finite(c.sigma * factor);
  params_.a = finite((sigma * sigma) / (1.0f / (float(c.bw) * float(c.bh))));
  params_.floor = (c.beta - 1) / c.beta;
  if (c.degrid != 0) {
    auto block = buffer<float>(fft.samples());
    const float peak = f.floating ? 1.0f : float((1 << f.bits) - 1);
    for (int y = 0; y < c.bh; ++y)
      for (int x = 0; x < c.bw; ++x)
        block[std::size_t(y) * c.bw + x] = (peak * wy_.analysis[y]) * wx_.analysis[x];
    grid_ = buffer<std::complex<float>>(fft.bins());
    fft.forward(block.data(), grid_.data());
    require(grid_[0].real() != 0, "FFT3D unusable grid DC");
  }
}
Plan::Plan(int w, int h, SampleFormat f, const DFTConfig& c)
    : geometry(geometry_dft(w, h, c)), format(f), algorithm(Algorithm::DFTTest), fft(c.block, c.block),
      kernel_(select_spectral(c.opt)), mean_scale_(c.zmean ? 1.0f : 0.0f), center_(c.mode == 0) {
  valid_format(f);
  auto win = dft_window(c.block, c.mode == 0 ? 0 : c.overlap, c.mode, c.swin, c.twin, c.sbeta, c.tbeta);
  h_ = std::move(win.h);
  const float scale = c.ftype < 2 ? win.wscale : 1.0f;
  params_ = {c.ftype,
             finite(c.sigma / scale),
             finite(c.sigma2 / scale),
             finite(c.pmin / win.wscale),
             finite(c.pmax / win.wscale),
             c.f0beta,
             0};
  if (c.zmean) {
    auto block = h_;
    for (float& v : block)
      v = finite(255.0f * v);
    grid_ = buffer<std::complex<float>>(fft.bins());
    fft.forward(block.data(), grid_.data());
    require(grid_[0].real() != 0, "DFTTest zmean requires nonzero template DC");
  }
}
template <class T>
void Plan::run(span2d::Plane<const T> src, span2d::Plane<T> dst) const {
  const auto& gx = geometry.x;
  const auto& gy = geometry.y;
  require(src.width() == gx.length && src.height() == gy.length && dst.width() == gx.length &&
              dst.height() == gy.length,
          "frame dimensions differ from plan");
  require(format.floating == std::is_same_v<T, float> && (format.floating || ((format.bits == 8) == (sizeof(T) == 1))),
          "frame storage differs from plan");
  const auto se = plane_extent<T>(src.width(), src.height(), src.stride_bytes());
  const auto de = plane_extent<T>(dst.width(), dst.height(), dst.stride_bytes());
  checked_plane(src.data(), src.width(), src.height(), src.stride_bytes(), se);
  checked_plane(dst.data(), dst.width(), dst.height(), dst.stride_bytes(), de);
  disjoint(src.data(), se, dst.data(), de);
  for (int y = 0; y < src.height(); ++y)
    for (int x = 0; x < src.width(); ++x)
      finite(float(src.row_ptr(y)[x]));
  auto accum = buffer<float>(mul_size(gx.cover, gy.cover));
  // One active block bounds FFT scratch independently of the number of origins.
  auto block = buffer<float>(fft.samples()), inverse = buffer<float>(fft.samples());
  auto spectrum = buffer<std::complex<float>>(fft.bins());
  auto row = algorithm == Algorithm::FFT3D ? buffer<float>(mul_size(gx.cover, gy.block)) : std::vector<float>{};
  const float base =
      !format.floating && format.chroma && algorithm == Algorithm::FFT3D ? float(1 << (format.bits - 1)) : 0;
  const float input_scale = format.floating ? 255.0f : 1.0f / float(1 << (format.bits - 8));
  const float volume = float(fft.width()) * float(fft.height());
  for (int by = 0; by < gy.count; ++by) {
    const int oy = by * gy.step;
    if (!row.empty())
      std::fill(row.begin(), row.end(), 0.0f);
    for (int bx = 0; bx < gx.count; ++bx) {
      const int ox = bx * gx.step;
      for (int y = 0; y < gy.block; ++y) {
        const auto* source = src.row_ptr(reflect(oy + y, gy));
        for (int x = 0; x < gx.block; ++x) {
          const auto i = std::size_t(y) * gx.block + x;
          const float q = float(source[reflect(ox + x, gx)]);
          block[i] = algorithm == Algorithm::FFT3D ? ((q - base) * wy_.analysis[y]) * wx_.analysis[x]
                                                   : finite(q * input_scale) * h_[i];
        }
      }
      fft.forward(block.data(), spectrum.data());
      const float scale = grid_.empty() ? 0 : finite(finite(mean_scale_ * spectrum[0].real()) / grid_[0].real());
      kernel_(spectrum.data(), grid_.empty() ? nullptr : grid_.data(), fft.bins(), scale, params_);
      fft.inverse(spectrum.data(), inverse.data());
      if (center_) {
        const int cy = gy.block / 2, cx = gx.block / 2;
        const auto i = std::size_t(cy) * gx.block + cx;
        accum[std::size_t(oy + cy) * gx.cover + ox + cx] = finite(finite(inverse[i] * volume) * h_[i]);
      } else
        for (int y = 0; y < gy.block; ++y)
          for (int x = 0; x < gx.block; ++x) {
            const auto i = std::size_t(y) * gx.block + x;
            if (algorithm == Algorithm::FFT3D) {
              auto& v = row[std::size_t(y) * gx.cover + ox + x];
              v = finite(v + finite(inverse[i] * wx_.synthesis[x]));
            } else {
              auto& v = accum[std::size_t(oy + y) * gx.cover + ox + x];
              v = finite(v + finite(finite(inverse[i] * volume) * h_[i]));
            }
          }
    }
    if (!row.empty())
      for (int y = 0; y < gy.block; ++y)
        for (int x = 0; x < gx.cover; ++x) {
          auto& v = accum[std::size_t(oy + y) * gx.cover + x];
          v = finite(v + finite(row[std::size_t(y) * gx.cover + x] * wy_.synthesis[y]));
        }
  }
  for (int y = 0; y < dst.height(); ++y)
    for (int x = 0; x < dst.width(); ++x) {
      const float z = accum[std::size_t(y + gy.offset) * gx.cover + x + gx.offset];
      if constexpr (std::is_same_v<T, float>) {
        dst.row_ptr(y)[x] = algorithm == Algorithm::FFT3D ? std::clamp(z, 0.0f, 1.0f) : finite(z * (1.0f / 255));
      } else {
        const float v = algorithm == Algorithm::FFT3D ? finite((z + 0.5f) + base)
                                                      : finite(finite(z * float(1 << (format.bits - 8))) + 0.5f);
        const float peak = float((1 << format.bits) - 1);
        dst.row_ptr(y)[x] = static_cast<T>(std::clamp(v, 0.0f, peak));
      }
    }
}
void Plan::process(span2d::Plane<const std::uint8_t> s, span2d::Plane<std::uint8_t> d) const {
  run(s, d);
}
void Plan::process(span2d::Plane<const std::uint16_t> s, span2d::Plane<std::uint16_t> d) const {
  run(s, d);
}
void Plan::process(span2d::Plane<const float> s, span2d::Plane<float> d) const {
  run(s, d);
}
} // namespace neo_fft
