#pragma once
#include "base/checked.hpp"
#include "spectral/fft_backend.hpp"
#include <complex>

namespace neo_fft {
// Distances and row strides are in elements. Each active slice is a Plane.
struct BatchLayout {
  std::size_t row_stride = 0, distance = 0, capacity = 0, active = 0;
};
class RealFFT {
public:
  RealFFT(int height, int width, FftProfile profile = FftProfile::native);
  int width() const { return width_; }
  int height() const { return height_; }
  int columns() const { return width_ / 2 + 1; }
  std::size_t samples() const { return std::size_t(width_) * height_; }
  std::size_t bins() const { return std::size_t(columns()) * height_; }
  FftProfile profile() const { return profile_; }
  int lanes() const { return backend_.lanes; }
  const char* backend_name() const { return backend_.name; }
  void forward(const float* in, BatchLayout real, std::complex<float>* out, BatchLayout spectrum) const;
  void inverse(const std::complex<float>* in, BatchLayout spectrum, float* out, BatchLayout real) const;
  void forward(const float* in, std::complex<float>* out) const;
  void inverse(const std::complex<float>* in, float* out) const;

private:
  int height_, width_;
  FftProfile profile_;
  const FftBackend& backend_;
};

class RealFFT3D {
public:
  RealFFT3D(int depth, int height, int width, FftProfile profile = FftProfile::native);
  int depth() const { return depth_; }
  int height() const { return height_; }
  int width() const { return width_; }
  int columns() const { return width_ / 2 + 1; }
  std::size_t samples() const { return std::size_t(depth_) * height_ * width_; }
  std::size_t bins() const { return std::size_t(depth_) * height_ * columns(); }
  FftProfile profile() const { return profile_; }
  int lanes() const { return backend_.lanes; }
  const char* backend_name() const { return backend_.name; }

  void forward(const float* in, std::complex<float>* out) const;
  void inverse(const std::complex<float>* in, float* out) const;

  void forward(const float* in, std::size_t batch, std::size_t in_dist,
               std::complex<float>* out, std::size_t out_dist) const;
  void inverse(const std::complex<float>* in, std::size_t batch, std::size_t in_dist,
               float* out, std::size_t out_dist) const;
  // Write only one H*W temporal plane per volume, starting at out + b*out_dist.
  // Invalid slices and unsupported shapes/profiles/group sizes return false without accessing buffers.
  bool try_inverse_slice(int slice, const std::complex<float>* in, std::size_t batch, std::size_t in_dist,
                          float* out, std::size_t out_dist) const;

private:
  int depth_, height_, width_;
  FftProfile profile_;
  const FftBackend& backend_;
};
} // namespace neo_fft
