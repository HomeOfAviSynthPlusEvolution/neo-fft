#pragma once
#include "algorithms/windows.hpp"
#include "algorithms/profile.hpp"
#include "algorithms/fft3d_profile.hpp"
#include <optional>
#include "kernels/spectral.hpp"
#include "kernels/spatial.hpp"
#include "runtime/workspace.hpp"
#include "runtime/workspace_pool.hpp"

namespace neo_fft {
enum class Algorithm { FFT3D, DFTTest };
struct SampleFormat {
  int bits = 8;
  bool floating = false, chroma = false;
};
struct FFT3DConfig {
  int bw = 32, bh = 32, ow = -1, oh = -1, wintype = 0, opt = 0, bt = 1;
  float sigma = 2, beta = 1, degrid = 1;
  std::optional<float> sigma2, sigma3, sigma4;
  EnhancementConfig enhancement;
};
struct DFTConfig {
  int block = 16, overlap = 12, mode = 1, swin = 0, twin = 7, ftype = 0, opt = 0, tbsize = 1;
  float sbeta = 2.5f, tbeta = 2.5f, sigma = 8, sigma2 = 8, pmin = 0, pmax = 500, f0beta = 1;
  bool zmean = true;
  DFTCurves curves;
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
  const int temporal_size = 1;
  const RealFFT fft;
  void process(span2d::Plane<const std::uint8_t> src, span2d::Plane<std::uint8_t> dst) const;
  void process(span2d::Plane<const std::uint16_t> src, span2d::Plane<std::uint16_t> dst) const;
  void process(span2d::Plane<const float> src, span2d::Plane<float> dst) const;

  void process(span2d::Plane<const std::uint8_t> src, span2d::Plane<std::uint8_t> dst, runtime::Workspace& ws) const;
  void process(span2d::Plane<const std::uint16_t> src, span2d::Plane<std::uint16_t> dst, runtime::Workspace& ws) const;
  void process(span2d::Plane<const float> src, span2d::Plane<float> dst, runtime::Workspace& ws) const;

  void process(span2d::Span<const span2d::Plane<const std::uint8_t>> sources, span2d::Plane<std::uint8_t> dst) const;
  void process(span2d::Span<const span2d::Plane<const std::uint16_t>> sources, span2d::Plane<std::uint16_t> dst) const;
  void process(span2d::Span<const span2d::Plane<const float>> sources, span2d::Plane<float> dst) const;

  void process(span2d::Span<const span2d::Plane<const std::uint8_t>> sources, span2d::Plane<std::uint8_t> dst,
               runtime::Workspace& ws) const;
  void process(span2d::Span<const span2d::Plane<const std::uint16_t>> sources, span2d::Plane<std::uint16_t> dst,
               runtime::Workspace& ws) const;
  void process(span2d::Span<const span2d::Plane<const float>> sources, span2d::Plane<float> dst,
               runtime::Workspace& ws) const;

  runtime::WorkspacePool& workspace_pool() const noexcept { return pool_; }

private:
  bool denoise_ = true;
  EnhancementTables enhancement_tables_;
  SpectralParams enhancement_params_;
  void enhance(std::complex<float>* spectrum) const;
  AxisWindow wx_, wy_;
  std::vector<float> h_;
  std::vector<float> h_synthesis_;
  std::vector<std::complex<float>> grid_;
  std::unique_ptr<RealFFT3D> fft3d_;
  std::vector<float> primary_;
  SpectralParams params_;
  SpectralKernel kernel_;
  Fft3dTemporalKernel temporal_kernel_ = nullptr;
  SpatialKernels spatial_;
  float mean_scale_ = 0;
  bool center_ = false;
  float sigma_eff_ = 0.0f;
  float norm_ = 0.0f;
  float beta_ = 1.0f;
  mutable runtime::WorkspacePool pool_;
  template <class T>
  void run(span2d::Span<const span2d::Plane<const T>> sources, span2d::Plane<T> dst, runtime::Workspace& ws) const;
};
} // namespace neo_fft
