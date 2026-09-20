#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/spectral.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/spectral.hpp"

HWY_BEFORE_NAMESPACE();
namespace neo_fft {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;
void Spectral(std::complex<float>* x, const std::complex<float>* grid, std::size_t count, float scale,
              const SpectralParams& p) {
  const hn::ScalableTag<float> d;
  if (!count)
    return;
  const auto lanes = hn::Lanes(d);
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0f), eps = hn::Set(d, 1e-15f);
  const auto a = hn::Set(d, p.a), b = hn::Set(d, p.b), low = hn::Set(d, p.low), high = hn::Set(d, p.high);
  auto check = [&](auto v) {
    if (!hn::AllTrue(d, hn::IsFinite(v)))
      throw std::runtime_error("non-finite sample or intermediate");
    return v;
  };
  // C++ explicitly permits float-array access to std::complex<float> storage.
  auto* output = reinterpret_cast<float*>(x);
  const auto* model = reinterpret_cast<const float*>(grid);
  std::size_t k = 0;
  for (; count - k >= lanes; k += lanes) {
    auto re = zero, im = zero, mr = zero, mi = zero;
    hn::LoadInterleaved2(d, output + 2 * k, re, im);
    if (grid) {
      hn::LoadInterleaved2(d, model + 2 * k, mr, mi);
      mr = check(hn::Mul(mr, hn::Set(d, scale)));
      mi = check(hn::Mul(mi, hn::Set(d, scale)));
    }
    re = check(hn::Sub(re, mr));
    im = check(hn::Sub(im, mi));
    const auto power = check(hn::Add(hn::Mul(re, re), hn::Mul(im, im)));
    auto gain = zero;
    switch (p.type) {
      case -1: {
        const auto q = hn::Add(power, eps);
        gain = hn::Max(hn::Div(hn::Sub(q, a), q), hn::Set(d, p.floor));
        break;
      }
      case 0: {
        gain = hn::Max(hn::Div(hn::Sub(power, a), hn::Add(power, eps)), zero);
        if (std::abs(p.exponent - 1.0f) < 0.00005f) {
        } else if (std::abs(p.exponent - .5f) < 0.00005f)
          gain = hn::Sqrt(gain);
        else {
          HWY_ALIGN float values[hn::MaxLanes(d)];
          hn::Store(gain, d, values);
          for (std::size_t i = 0; i < lanes; ++i)
            values[i] = std::pow(values[i], p.exponent);
          gain = hn::Load(d, values);
        }
        break;
      }
      case 1:
        gain = hn::IfThenElse(hn::Lt(power, a), zero, one);
        break;
      case 2:
        gain = a;
        break;
      case 3:
        gain = hn::IfThenElse(hn::And(hn::Ge(power, low), hn::Le(power, high)), a, b);
        break;
      case 4: {
        const auto q = hn::Add(power, eps);
        const auto num = check(hn::Mul(q, high));
        const auto den = check(hn::Mul(check(hn::Add(q, low)), check(hn::Add(q, high))));
        gain = hn::Mul(a, hn::Sqrt(check(hn::Div(num, den))));
        break;
      }
      default:
        throw std::invalid_argument("invalid spectral filter type");
    }
    check(gain);
    re = check(hn::Add(hn::Mul(gain, re), mr));
    im = check(hn::Add(hn::Mul(gain, im), mi));
    hn::StoreInterleaved2(re, im, d, output + 2 * k);
  }
  spectral_scalar(x + k, grid ? grid + k : nullptr, count - k, scale, p);
}
const char* Target() {
  return hwy::TargetName(HWY_TARGET);
}
} // namespace HWY_NAMESPACE
} // namespace neo_fft
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(Spectral);
HWY_EXPORT(Target);
SpectralKernel select_spectral(int opt) {
  require(opt == 0 || opt == 1, "unsupported opt: expected 0 or 1");
  return opt == 1 ? spectral_scalar : HWY_DYNAMIC_DISPATCH(Spectral);
}
const char* spectral_target(int opt) {
  select_spectral(opt);
  return opt == 1 ? "scalar" : HWY_DYNAMIC_DISPATCH(Target)();
}
} // namespace neo_fft
#endif
