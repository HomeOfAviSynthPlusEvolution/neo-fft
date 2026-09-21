#pragma once
#include "algorithms/windows.hpp"
#include "kernels/spectral.hpp"
#include "runtime/workspace.hpp"
#include "runtime/workspace_pool.hpp"

namespace neo_fft {
enum class Algorithm { FFT3D, DFTTest };
struct SampleFormat {
  int bits = 8;
  bool floating = false, chroma = false;
};
struct FFT3DConfig {
  int bw = 32, bh = 32, ow = -1, oh = -1, wintype = 0, opt = 0;
  float sigma = 2, beta = 1, degrid = 1;
};
struct DFTConfig {
  int block = 16, overlap = 12, mode = 1, swin = 0, twin = 7, ftype = 0, opt = 0;
  float sbeta = 2.5f, tbeta = 2.5f, sigma = 8, sigma2 = 8, pmin = 0, pmax = 500, f0beta = 1;
  bool zmean = true;
};
void validate(const FFT3DConfig& c);
void validate(const DFTConfig& c);
class Plan {
public:
  Plan(int width, int height, SampleFormat format, const FFT3DConfig& config);
  Plan(int width, int height, SampleFormat format, const DFTConfig& config);
  const Geometry geometry;
  const SampleFormat format;
  const Algorithm algorithm;
  const RealFFT fft;
  void process(span2d::Plane<const std::uint8_t> src, span2d::Plane<std::uint8_t> dst) const;
  void process(span2d::Plane<const std::uint16_t> src, span2d::Plane<std::uint16_t> dst) const;
  void process(span2d::Plane<const float> src, span2d::Plane<float> dst) const;

  void process(span2d::Plane<const std::uint8_t> src, span2d::Plane<std::uint8_t> dst, runtime::Workspace& ws) const;
  void process(span2d::Plane<const std::uint16_t> src, span2d::Plane<std::uint16_t> dst, runtime::Workspace& ws) const;
  void process(span2d::Plane<const float> src, span2d::Plane<float> dst, runtime::Workspace& ws) const;

  runtime::WorkspacePool& workspace_pool() const noexcept { return pool_; }

private:
  AxisWindow wx_, wy_;
  std::vector<float> h_;
  std::vector<std::complex<float>> grid_;
  SpectralParams params_;
  SpectralKernel kernel_;
  float mean_scale_ = 0;
  bool center_ = false;
  mutable runtime::WorkspacePool pool_;
  template <class T>
  void run(span2d::Plane<const T> src, span2d::Plane<T> dst, runtime::Workspace& ws) const;
};
} // namespace neo_fft
