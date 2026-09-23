#pragma once
#include "algorithms/windows.hpp"
#include "algorithms/temporal_grid.hpp"
#include "algorithms/profile.hpp"
#include "algorithms/dft_noise.hpp"
#include "algorithms/fft3d_profile.hpp"
#include <optional>
#include "runtime/published_model.hpp"
#include "kernels/spectral.hpp"
#include "kernels/spatial.hpp"
#include "kernels/model.hpp"
#include "kernels/kalman.hpp"
#include "kernels/rows.hpp"
#include "runtime/executor.hpp"
#include "runtime/workspace.hpp"
#include "runtime/workspace_pool.hpp"
#include "runtime/spectra_cache.hpp"

namespace neo_fft {
enum class Algorithm { FFT3D, DFTTest };
struct SampleFormat {
  int bits = 8;
  bool floating = false, chroma = false;
};
struct FFT3DConfig {
  int bw = 32, bh = 32, ow = -1, oh = -1, wintype = 0, opt = 0, bt = 1;
  float sigma = 2, beta = 1, degrid = 1, kratio = 2;
  std::optional<float> sigma2, sigma3, sigma4;
  EnhancementConfig enhancement;
  float pfactor = 0, pcutoff = .1f;
  int pframe = 0, px = 0, py = 0;
  bool pshow = false;
  int left = 0, top = 0, right = 0, bottom = 0;
  bool interlaced = false, mt = false;
  int ncpu = 2;
  int cache_frames = -1, cache_mb = runtime::SpectraCache::default_mb;
};
struct DFTConfig {
  int block = 16, overlap = 12, mode = 1, swin = 0, twin = 7, ftype = 0, opt = 0, tbsize = 1;
  float sbeta = 2.5f, tbeta = 2.5f, sigma = 8, sigma2 = 8, pmin = 0, pmax = 500, f0beta = 1;
  bool zmean = true;
  int temporal_mode = 0, temporal_overlap = 0;
  int dither = 0, dither_seed = 0, threads = 1, fft_threads = 1;
  DFTCurves curves;
  std::vector<NoiseLocation> locations;
  std::optional<float> alpha;
};
void validate(const FFT3DConfig& c);
void validate(const DFTConfig& c);
class Plan {
public:
  Plan(int width, int height, SampleFormat format, const FFT3DConfig& config, std::shared_ptr<runtime::Executor> executor = {}, std::shared_ptr<runtime::Retention> retention = {});
  Plan(int width, int height, SampleFormat format, const DFTConfig& config, std::shared_ptr<const DFTNoise> noise = {}, std::shared_ptr<runtime::Executor> executor = {}, std::shared_ptr<runtime::Retention> retention = {});
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

  template<class T> void process_at(span2d::Span<const span2d::Plane<const T>> sources,
                                    span2d::Plane<T> dst,int frame,int plane,span2d::Span<const int> targets = {},runtime::SpectraCache* cache = nullptr,runtime::SpectraCache::Request* registration=nullptr) const {
    auto lease=pool_.acquire(); run(sources,dst,*lease,frame,plane,nullptr,targets,cache,registration);
  }
  std::shared_ptr<const DFTNoise> dft_noise() const { return dft_noise_; }
  const runtime::Executor& executor() const {return *executor_;}
  CopyRow copy_row() const {return copy_row_;}
  bool kalman() const { return kalman_; }
  KalmanState initial_kalman() const;
  template<class T> void advance_kalman(span2d::Plane<const T> source, KalmanState& state) const;
  template<class T> void render_kalman(span2d::Plane<const T> source,span2d::Plane<T> dst,const KalmanState& state) const {
    auto lease=pool_.acquire(); run<T>({&source,1},dst,*lease,0,0,&state);
  }
  bool preview() const { return preview_; }
  bool needs_pattern_frame() const { return sampled_ && denoise_; }
  bool pattern_ready() const { return bool(sampled_model_.get()); }
  runtime::PublishedModel::Model pattern_power() const { return sampled_model_.get(); }
  template<class T> runtime::PublishedModel::Model pattern_candidate(span2d::Plane<const T> source) const;
  void publish_pattern(runtime::PublishedModel::Model model) const { sampled_model_.publish(std::move(model)); }
  void prepare_pattern(span2d::Plane<const std::uint8_t> source) const;
  void prepare_pattern(span2d::Plane<const std::uint16_t> source) const;
  void prepare_pattern(span2d::Plane<const float> source) const;
  runtime::WorkspacePool& workspace_pool() const noexcept { return pool_; }

private:
  std::shared_ptr<runtime::Executor> executor_;
  KalmanKernel kalman_kernel_ = kalman_scalar;
  CopyRow copy_row_ = copy_row_scalar;
  DitherNoise dither_noise_ = dither_noise_scalar;
  bool kalman_ = false;
  float kalman_r0_ = 0, kalman_ratio2_ = 0;
  std::size_t state_bins() const;
  bool temporal_ola_ = false;
  int dither_ = 0, dither_seed_ = 0;
  std::shared_ptr<const DFTNoise> dft_noise_;
  bool denoise_ = true, sampled_ = false, preview_ = false;
  int px_ = 0, py_ = 0;
  float pfactor_ = 0;
  std::vector<float> sample_weights_;
  mutable runtime::PublishedModel sampled_model_;
  template<class T> std::pair<int,int> select_pattern(span2d::Plane<const T> source, runtime::Workspace& ws, std::vector<float>* power) const;
  template<class T> void build_pattern(span2d::Plane<const T> source) const;
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
  ModelKernels model_;
  template<class T> void gather_pattern(span2d::Plane<const T>,float*,int,int,bool) const;
  float mean_scale_ = 0;
  bool center_ = false;
  float sigma_eff_ = 0.0f;
  float norm_ = 0.0f;
  float beta_ = 1.0f;
  mutable runtime::WorkspacePool pool_;
  template <class T>
  void run(span2d::Span<const span2d::Plane<const T>> sources, span2d::Plane<T> dst, runtime::Workspace& ws, int frame=0, int plane=0, const KalmanState* kalman=nullptr,span2d::Span<const int> targets={},runtime::SpectraCache* cache=nullptr,runtime::SpectraCache::Request* registration=nullptr) const;
};
} // namespace neo_fft
