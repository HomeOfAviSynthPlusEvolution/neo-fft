#pragma once
#include "spectral/fft.hpp"

namespace neo_fft {
enum class PrimaryMode { Uniform, Table };
struct SpectralParams {
  // -1: FFT3D; 0..4: DFTTest. All constants are already calibrated.
  int type = -1;
  float a = 0, b = 0, low = 0, high = 0, exponent = 1, floor = 0;
  PrimaryMode primary_mode = PrimaryMode::Uniform;
  // Immutable calibrated logical bins, disjoint from writable spectrum.
  span2d::Span<const float> primary{};
};
using SpectralKernel = void (*)(std::complex<float>*, const std::complex<float>*, std::size_t, float,
                                const SpectralParams&);
void spectral_scalar(std::complex<float>* spectrum, const std::complex<float>* grid, std::size_t count,
                     float grid_scale, const SpectralParams& p);
struct Fft3dTwiddles {
  std::complex<float> fwd[6][5][5]{};
  std::complex<float> inv[6][5][5]{};
};
const Fft3dTwiddles& get_fft3d_twiddles() noexcept;

using Fft3dTemporalKernel = void (*)(const std::complex<float>* const* spectra, int T, int c, std::size_t bins,
                                     float degrid, const std::complex<float>* grid, float noise, float lower,
                                     std::complex<float>* out);
void fft3d_temporal_scalar(const std::complex<float>* const* spectra, int T, int c, std::size_t bins,
                           float degrid, const std::complex<float>* grid, float noise, float lower,
                           std::complex<float>* out);
Fft3dTemporalKernel select_fft3d_temporal(int opt);
const char* fft3d_temporal_target(int opt);

void fft3d_temporal_filter(const std::complex<float>* const* spectra, int T, int c, std::size_t bins,
                           float degrid, const std::complex<float>* grid, float noise, float lower,
                           std::complex<float>* out);

SpectralKernel select_spectral(int opt);
const char* spectral_target(int opt);
std::size_t optimal_l2_working_set_bytes() noexcept;
int optimal_simd_lanes() noexcept;
} // namespace neo_fft
