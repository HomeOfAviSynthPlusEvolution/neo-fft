#pragma once
#include "spectral/fft.hpp"
#include "kernels/model.hpp"
#include "runtime/published_model.hpp"
#include <functional>

namespace neo_fft {
struct DFTConfig;
struct NoiseLocation { int frame, plane, y, x; };
// Shared by all output planes of one filter. The host fills one chronological
// spatial slice in 8-bit amplitude units; no host frame survives that callback.
class DFTNoise {
public:
  explicit DFTNoise(const DFTConfig& config);
  using Gather = std::function<void(const NoiseLocation&, int, span2d::Span<float>)>;
  void prepare(const Gather& gather) const;
  const ModelKernels& kernels() const { return kernels_; }
  runtime::PublishedModel::Model power() const { return model_.get(); }
  const std::vector<NoiseLocation> locations;
  const int temporal_size, block_size;
  std::size_t working_set_bytes() const { return working_set_bytes_; }
private:
  RealFFT3D fft_;
  ModelKernels kernels_;
  std::size_t working_set_bytes_;
  std::vector<float> window_;
  std::vector<std::complex<float>> grid_;
  float calibration_;
  bool zmean_;
  mutable runtime::PublishedModel model_;
};
} // namespace neo_fft
