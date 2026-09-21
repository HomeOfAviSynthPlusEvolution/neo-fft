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
  auto compute_gain = [&](auto power) {
    switch (p.type) {
      case -1: {
        const auto q = hn::Add(power, eps);
        return hn::Max(hn::Div(hn::Sub(q, a), q), hn::Set(d, p.floor));
      }
      case 0: {
        auto g = hn::Max(hn::Div(hn::Sub(power, a), hn::Add(power, eps)), zero);
        if (std::abs(p.exponent - 1.0f) < 0.00005f) {
          return g;
        } else if (std::abs(p.exponent - .5f) < 0.00005f) {
          return hn::Sqrt(g);
        } else {
          HWY_ALIGN float values[hn::MaxLanes(d)];
          hn::Store(g, d, values);
          for (std::size_t i = 0; i < lanes; ++i)
            values[i] = std::pow(values[i], p.exponent);
          return hn::Load(d, values);
        }
      }
      case 1:
        return hn::IfThenElse(hn::Lt(power, a), zero, one);
      case 2:
        return a;
      case 3:
        return hn::IfThenElse(hn::And(hn::Ge(power, low), hn::Le(power, high)), a, b);
      case 4: {
        const auto q = hn::Add(power, eps);
        const auto num = hn::Mul(q, high);
        const auto den = hn::Mul(hn::Add(q, low), hn::Add(q, high));
        return hn::Mul(a, hn::Sqrt(hn::Div(num, den)));
      }
      default:
        throw std::invalid_argument("invalid spectral filter type");
    }
  };

  // C++ explicitly permits float-array access to std::complex<float> storage.
  auto* output = reinterpret_cast<float*>(x);
  const auto* model = reinterpret_cast<const float*>(grid);
  auto non_finite = hn::MaskFalse(d);
  std::size_t k = 0;
  if (!grid) {
    for (; count - k >= lanes; k += lanes) {
      auto re = zero, im = zero;
      hn::LoadInterleaved2(d, output + 2 * k, re, im);
      const auto power = hn::Add(hn::Mul(re, re), hn::Mul(im, im));
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(power)));
      const auto gain = compute_gain(power);
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(gain)));
      re = hn::Mul(gain, re);
      im = hn::Mul(gain, im);
      hn::StoreInterleaved2(re, im, d, output + 2 * k);
    }
  } else {
    const auto scale_v = hn::Set(d, scale);
    for (; count - k >= lanes; k += lanes) {
      auto re = zero, im = zero, mr = zero, mi = zero;
      hn::LoadInterleaved2(d, output + 2 * k, re, im);
      hn::LoadInterleaved2(d, model + 2 * k, mr, mi);
      mr = hn::Mul(mr, scale_v);
      mi = hn::Mul(mi, scale_v);
      re = hn::Sub(re, mr);
      im = hn::Sub(im, mi);
      const auto power = hn::Add(hn::Mul(re, re), hn::Mul(im, im));
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(power)));
      const auto gain = compute_gain(power);
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(gain)));
      re = hn::Add(hn::Mul(gain, re), mr);
      im = hn::Add(hn::Mul(gain, im), mi);
      hn::StoreInterleaved2(re, im, d, output + 2 * k);
    }
  }
  if (!hn::AllFalse(d, non_finite))
    throw std::runtime_error("non-finite sample or intermediate");
  spectral_scalar(x + k, grid ? grid + k : nullptr, count - k, scale, p);
}
const char* Target() {
  return hwy::TargetName(HWY_TARGET);
}
int SimdLanes() {
  const hn::ScalableTag<float> d;
  return static_cast<int>(hn::Lanes(d));
}
} // namespace HWY_NAMESPACE
} // namespace neo_fft
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(Spectral);
HWY_EXPORT(Target);
HWY_EXPORT(SimdLanes);
SpectralKernel select_spectral(int opt) {
  require(opt == 0 || opt == 1, "unsupported opt: expected 0 or 1");
  return opt == 1 ? spectral_scalar : HWY_DYNAMIC_DISPATCH(Spectral);
}
const char* spectral_target(int opt) {
  select_spectral(opt);
  return opt == 1 ? "scalar" : HWY_DYNAMIC_DISPATCH(Target)();
}
int optimal_simd_lanes() noexcept {
  return HWY_DYNAMIC_DISPATCH(SimdLanes)();
}
std::size_t optimal_l2_working_set_bytes() noexcept {
  const int64_t targets = hwy::SupportedTargets();
#if defined(HWY_AVX3_ZEN4) || defined(HWY_AVX3_SPR)
  constexpr int64_t kBigL2 = (HWY_AVX3_ZEN4 | HWY_AVX3_SPR);
  if (targets & kBigL2) {
    return 512 * 1024;
  }
#endif
#if defined(HWY_AVX3) || defined(HWY_AVX3_DL)
  constexpr int64_t kMidL2 = (HWY_AVX3 | HWY_AVX3_DL);
  if (targets & kMidL2) {
    return 384 * 1024;
  }
#endif
#if defined(HWY_AVX2)
  if (targets & HWY_AVX2) {
    return 256 * 1024;
  }
#endif
  return 128 * 1024;
}
} // namespace neo_fft
#endif
