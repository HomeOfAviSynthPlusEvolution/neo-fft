#pragma once
#include "base/checked.hpp"
#include <complex>
namespace neo_fft {
struct KalmanState {
  std::vector<std::complex<float>> last;
  // Canonical states have equal real/imaginary covariance components. C and Q
  // remain independent; each stores the shared component once, without loss.
  std::vector<float> covariance, process;
  std::size_t bytes() const {
    return add_size(sizeof(*this),add_size(mul_size(last.capacity(),sizeof(std::complex<float>)),
        mul_size(add_size(covariance.capacity(),process.capacity()),sizeof(float))));
  }
};
// Update last/covariance/process in place. Their count-element ranges must be
// mutually disjoint and must not overlap x or pattern. On failure the private
// request state may be partially updated and must be discarded.
using KalmanKernel = void (*)(const std::complex<float>*, std::complex<float>*, float*,
                            float*, const float*, float, float, std::size_t);
KalmanKernel select_kalman(int opt);
void kalman_scalar(const std::complex<float>* x, std::complex<float>* last, float* covariance,
                   float* process, const float* pattern, float uniform, float ratio2, std::size_t count);
} // namespace neo_fft
