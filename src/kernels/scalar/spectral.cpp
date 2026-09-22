#include "kernels/spectral.hpp"
#include <algorithm>

namespace neo_fft {
void spectral_scalar(std::complex<float>* x, const std::complex<float>* grid, std::size_t count, float scale,
                     const SpectralParams& p) {
  require(p.primary_mode != PrimaryMode::Table || p.primary.size() == count, "primary table shape mismatch");
  require(!p.enhancement.sharpen || p.enhancement.sharpen_window.size() == count, "sharpen table shape mismatch");
  require(!p.enhancement.dehalo || p.enhancement.halo_window.size() == count, "halo table shape mismatch");
  for (std::size_t k = 0; k < count; ++k) {
    const float a = p.primary_mode == PrimaryMode::Table ? p.primary[k] : p.a;
    const float mr = grid ? finite(scale * grid[k].real()) : 0;
    const float mi = grid ? finite(scale * grid[k].imag()) : 0;
    const float re = finite(x[k].real() - mr), im = finite(x[k].imag() - mi);
    const float power = finite(re * re + im * im);
    float gain = 0;
    switch (p.type) {
      case -2:
        gain = 1;
        break;
      case -1: {
        const float q = power + 1e-15f;
        gain = std::max((q - a) / q, p.floor);
        break;
      }
      case 0: {
        const float wiener = std::max((power - a) / (power + 1e-15f), 0.0f);
        gain = std::abs(p.exponent - 1.0f) < 0.00005f   ? wiener
               : std::abs(p.exponent - 0.5f) < 0.00005f ? std::sqrt(wiener)
                                                        : std::pow(wiener, p.exponent);
        break;
      }
      case 1:
        gain = power < a ? 0.0f : 1.0f;
        break;
      case 2:
        gain = a;
        break;
      case 3:
        gain = power >= p.low && power <= p.high ? a : p.b;
        break;
      case 4: {
        const float q = power + 1e-15f;
        const float num = finite(q * p.high), den = finite(finite(q + p.low) * finite(q + p.high));
        gain = a * std::sqrt(finite(num / den));
        break;
      }
      default:
        throw std::invalid_argument("invalid spectral filter type");
    }
    if (p.enhancement.active())
      gain = finite(gain * enhancement_gain(finite(power + 1e-15f), k, p.enhancement));
    finite(gain);
    x[k] = {finite(gain * re + mr), finite(gain * im + mi)};
  }
}

const Fft3dTwiddles& get_fft3d_twiddles() noexcept {
  static const Fft3dTwiddles twiddles = []() {
    Fft3dTwiddles t{};
    constexpr double kPi = 3.1415926535897932384626433832795;
    for (int T = 2; T <= 5; ++T) {
      for (int m = 0; m < T; ++m) {
        for (int j = 0; j < T; ++j) {
          const double angle = -2.0 * kPi * double(j * m) / double(T);
          t.fwd[T][m][j] = {float(std::cos(angle)), float(std::sin(angle))};
        }
      }
      for (int c = 0; c < T; ++c) {
        for (int m = 0; m < T; ++m) {
          const double inv_angle = +2.0 * kPi * double(c * m) / double(T);
          t.inv[T][c][m] = {float(std::cos(inv_angle)), float(std::sin(inv_angle))};
        }
      }
    }
    return t;
  }();
  return twiddles;
}

void fft3d_temporal_scalar(const std::complex<float>* const* spectra, int T, int c, std::size_t bins,
                           float degrid, const std::complex<float>* grid, NoisePower noise, float lower,
                           std::complex<float>* out) {
  require(noise.mode != PrimaryMode::Table || noise.table.size() == bins, "noise table shape mismatch");
  require(T >= 1 && T <= 5, "FFT3D invalid T");
  require(c >= 0 && c < T, "FFT3D invalid c");
  if (T == 1) {
    const float mr_scale = (grid && degrid != 0.0f && grid[0].real() != 0.0f)
                               ? finite((degrid * spectra[0][0].real()) / grid[0].real())
                               : 0.0f;
    for (std::size_t k = 0; k < bins; ++k) {
      const float mr = grid ? finite(mr_scale * grid[k].real()) : 0.0f;
      const float mi = grid ? finite(mr_scale * grid[k].imag()) : 0.0f;
      const float re = finite(spectra[0][k].real() - mr);
      const float im = finite(spectra[0][k].imag() - mi);
      const float power = finite(re * re + im * im);
      const float q = power + 1e-15f;
      const float gain = std::max((q - noise.at(k)) / q, lower);
      out[k] = {finite(gain * re + mr), finite(gain * im + mi)};
    }
    return;
  }

  const auto& twiddles = get_fft3d_twiddles();
  const auto* fwd_twiddle = twiddles.fwd[T];
  const auto* inv_twiddle = twiddles.inv[T][c];

  const float g_ratio = (grid && degrid != 0.0f && grid[0].real() != 0.0f)
                            ? finite((degrid * spectra[c][0].real()) / grid[0].real())
                            : 0.0f;

  const float inv_T = 1.0f / float(T);

  for (std::size_t k = 0; k < bins; ++k) {
    const std::complex<float> M = grid ? (g_ratio * grid[k]) : std::complex<float>(0.0f, 0.0f);
    const std::complex<float> gridT = M * float(T);

    std::complex<float> F[5]{};
    for (int m = 0; m < T; ++m) {
      std::complex<float> sum(0.0f, 0.0f);
      for (int j = 0; j < T; ++j) {
        sum += spectra[j][k] * fwd_twiddle[m][j];
      }
      F[m] = sum;
    }

    std::complex<float> R[5];
    R[0] = F[0] - gridT;
    for (int m = 1; m < T; ++m) {
      R[m] = F[m];
    }

    std::complex<float> R_filtered[5];
    for (int m = 0; m < T; ++m) {
      const float power = finite(R[m].real() * R[m].real() + R[m].imag() * R[m].imag());
      const float q = power + 1e-15f;
      const float gain = std::max((q - noise.at(k)) / q, lower);
      R_filtered[m] = gain * R[m];
    }

    std::complex<float> F_out[5];
    F_out[0] = R_filtered[0] + gridT;
    for (int m = 1; m < T; ++m) {
      F_out[m] = R_filtered[m];
    }

    std::complex<float> y(0.0f, 0.0f);
    for (int m = 0; m < T; ++m) {
      y += F_out[m] * inv_twiddle[m];
    }
    out[k] = y * inv_T;
  }
}

void fft3d_temporal_filter(const std::complex<float>* const* spectra, int T, int c, std::size_t bins,
                           float degrid, const std::complex<float>* grid, NoisePower noise, float lower,
                           std::complex<float>* out) {
  select_fft3d_temporal(0)(spectra, T, c, bins, degrid, grid, noise, lower, out);
}

#if !NEO_FFT_ENABLE_HIGHWAY
SpectralKernel select_spectral(int opt) {
  return spectral_scalar;
}
const char* spectral_target(int opt) {
  select_spectral(opt);
  return "scalar (Highway disabled)";
}
Fft3dTemporalKernel select_fft3d_temporal(int opt) {
  return fft3d_temporal_scalar;
}
const char* fft3d_temporal_target(int opt) {
  select_fft3d_temporal(opt);
  return "scalar (Highway disabled)";
}
std::size_t optimal_l2_working_set_bytes() noexcept {
  return 256 * 1024;
}
int optimal_simd_lanes() noexcept {
  return 4;
}
#endif
} // namespace neo_fft
