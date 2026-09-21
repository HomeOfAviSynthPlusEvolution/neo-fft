#include "kernels/spectral.hpp"
#include <algorithm>

namespace neo_fft {
void spectral_scalar(std::complex<float>* x, const std::complex<float>* grid, std::size_t count, float scale,
                     const SpectralParams& p) {
  for (std::size_t k = 0; k < count; ++k) {
    const float mr = grid ? finite(scale * grid[k].real()) : 0;
    const float mi = grid ? finite(scale * grid[k].imag()) : 0;
    const float re = finite(x[k].real() - mr), im = finite(x[k].imag() - mi);
    const float power = finite(re * re + im * im);
    float gain = 0;
    switch (p.type) {
      case -1: {
        const float q = power + 1e-15f;
        gain = std::max((q - p.a) / q, p.floor);
        break;
      }
      case 0: {
        const float a = std::max((power - p.a) / (power + 1e-15f), 0.0f);
        gain = std::abs(p.exponent - 1.0f) < 0.00005f   ? a
               : std::abs(p.exponent - 0.5f) < 0.00005f ? std::sqrt(a)
                                                        : std::pow(a, p.exponent);
        break;
      }
      case 1:
        gain = power < p.a ? 0.0f : 1.0f;
        break;
      case 2:
        gain = p.a;
        break;
      case 3:
        gain = power >= p.low && power <= p.high ? p.a : p.b;
        break;
      case 4: {
        const float q = power + 1e-15f;
        const float num = finite(q * p.high), den = finite(finite(q + p.low) * finite(q + p.high));
        gain = p.a * std::sqrt(finite(num / den));
        break;
      }
      default:
        throw std::invalid_argument("invalid spectral filter type");
    }
    finite(gain);
    x[k] = {finite(gain * re + mr), finite(gain * im + mi)};
  }
}
#if !NEO_FFT_ENABLE_HIGHWAY
SpectralKernel select_spectral(int opt) {
  require(opt == 0 || opt == 1, "unsupported opt: expected 0 or 1");
  return spectral_scalar;
}
const char* spectral_target(int opt) {
  select_spectral(opt);
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
