#pragma once
#include "spectral/fft.hpp"

namespace neo_fft {
struct SpectralParams {
  // -1: FFT3D; 0..4: DFTTest. All constants are already calibrated.
  int type = -1;
  float a = 0, b = 0, low = 0, high = 0, exponent = 1, floor = 0;
};
using SpectralKernel = void (*)(std::complex<float>*, const std::complex<float>*, std::size_t, float,
                                const SpectralParams&);
void spectral_scalar(std::complex<float>* spectrum, const std::complex<float>* grid, std::size_t count,
                     float grid_scale, const SpectralParams& p);
SpectralKernel select_spectral(int opt);
const char* spectral_target(int opt);
std::size_t optimal_l2_working_set_bytes() noexcept;
int optimal_simd_lanes() noexcept;
} // namespace neo_fft
