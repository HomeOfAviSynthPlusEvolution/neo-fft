#pragma once
#include "base/checked.hpp"
#include <complex>
namespace neo_fft {
struct KalmanState {
  std::vector<std::complex<float>> last, covariance, process;
  std::size_t bytes() const {
    return add_size(sizeof(*this),mul_size(add_size(add_size(last.capacity(),covariance.capacity()),process.capacity()),sizeof(std::complex<float>)));
  }
};
using KalmanKernel = void (*)(const std::complex<float>*, std::complex<float>*, std::complex<float>*,
                            std::complex<float>*, const float*, float, float, std::size_t);
void kalman_scalar(const std::complex<float>* x, std::complex<float>* last, std::complex<float>* covariance,
                   std::complex<float>* process, const float* pattern, float uniform, float ratio2, std::size_t count);
} // namespace neo_fft
