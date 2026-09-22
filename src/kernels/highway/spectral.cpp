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
  require(p.primary_mode != PrimaryMode::Table || p.primary.size() == count, "primary table shape mismatch");
  require(!p.enhancement.sharpen || p.enhancement.sharpen_window.size() == count, "sharpen table shape mismatch");
  require(!p.enhancement.dehalo || p.enhancement.halo_window.size() == count, "halo table shape mismatch");
  const hn::ScalableTag<float> d;
  if (!count)
    return;
  const auto lanes = hn::Lanes(d);
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0f), eps = hn::Set(d, 1e-15f);
  const auto uniform_a = hn::Set(d, p.a), b = hn::Set(d, p.b), low = hn::Set(d, p.low), high = hn::Set(d, p.high);
  auto non_finite = hn::MaskFalse(d);
  auto compute_gain = [&](auto power, std::size_t offset) {
    const auto a = p.primary_mode == PrimaryMode::Table ? hn::LoadU(d, p.primary.data() + offset) : uniform_a;
    switch (p.type) {
      case -2:
        return one;
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
          constexpr hn::ScalableTag<float> df;
          HWY_ALIGN float values[hn::MaxLanes(df)];
          hn::Store(g, df, values);
          for (std::size_t i = 0; i < lanes; ++i)
            values[i] = std::pow(values[i], p.exponent);
          return hn::Load(df, values);
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
        non_finite = hn::Or(non_finite, hn::Not(hn::And(hn::IsFinite(num), hn::IsFinite(den))));
        return hn::Mul(a, hn::Sqrt(hn::Div(num, den)));
      }
      default:
        throw std::invalid_argument("invalid spectral filter type");
    }
  };

  // C++ explicitly permits float-array access to std::complex<float> storage.
  auto* output = reinterpret_cast<float*>(x);
  const auto* model = reinterpret_cast<const float*>(grid);
  auto enhance = [&](auto power, std::size_t offset) {
    const auto& e = p.enhancement;
    auto check = [&](auto v) { non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(v))); return v; };
    auto gain = one;
    const auto q = check(hn::Add(power, eps));
    if (e.sharpen != 0 && e.b != 0) {
      const auto num = check(hn::Mul(q, hn::Set(d, e.b)));
      const auto den = check(hn::Mul(check(hn::Add(q, hn::Set(d, e.a))), check(hn::Add(q, hn::Set(d, e.b)))));
      const auto strength = check(hn::Mul(hn::Set(d, e.sharpen), hn::LoadU(d, e.sharpen_window.data() + offset)));
      gain = check(hn::Add(one, check(hn::Mul(strength, check(hn::Sqrt(check(hn::Div(num, den))))))));
    }
    if (e.dehalo != 0) {
      const auto qc = check(hn::Add(q, hn::Set(d, e.c)));
      const auto suppression = check(hn::Mul(check(hn::Mul(hn::Set(d, e.dehalo), hn::LoadU(d, e.halo_window.data() + offset))), q));
      gain = check(hn::Mul(gain, check(hn::Div(qc, check(hn::Add(qc, suppression))))));
    }
    return gain;
  };
  std::size_t k = 0;
  if (!grid) {
    for (; count - k >= lanes; k += lanes) {
      auto re = zero, im = zero;
      hn::LoadInterleaved2(d, output + 2 * k, re, im);
      const auto power = hn::Add(hn::Mul(re, re), hn::Mul(im, im));
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(power)));
      auto gain = compute_gain(power, k);
      if (p.enhancement.active()) gain = hn::Mul(gain, enhance(power, k));
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(gain)));
      re = hn::Mul(gain, re);
      im = hn::Mul(gain, im);
      non_finite = hn::Or(non_finite, hn::Not(hn::And(hn::IsFinite(re), hn::IsFinite(im))));
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
      auto gain = compute_gain(power, k);
      if (p.enhancement.active()) gain = hn::Mul(gain, enhance(power, k));
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(gain)));
      re = hn::Add(hn::Mul(gain, re), mr);
      im = hn::Add(hn::Mul(gain, im), mi);
      non_finite = hn::Or(non_finite, hn::Not(hn::And(hn::IsFinite(re), hn::IsFinite(im))));
      hn::StoreInterleaved2(re, im, d, output + 2 * k);
    }
  }
  if (!hn::AllFalse(d, non_finite))
    throw std::runtime_error("non-finite sample or intermediate");
  auto tail = p;
  if (p.primary_mode == PrimaryMode::Table)
    tail.primary = {p.primary.data() + k, count - k};
  if (p.enhancement.sharpen) tail.enhancement.sharpen_window = {p.enhancement.sharpen_window.data() + k, count - k};
  if (p.enhancement.dehalo) tail.enhancement.halo_window = {p.enhancement.halo_window.data() + k, count - k};
  spectral_scalar(x + k, grid ? grid + k : nullptr, count - k, scale, tail);
}

void Fft3dTemporal(const std::complex<float>* const* spectra, int T, int c, std::size_t bins,
                   float degrid, const std::complex<float>* grid, NoisePower noise, float lower,
                   std::complex<float>* out) {
  require(noise.mode != PrimaryMode::Table || noise.table.size() == bins, "noise table shape mismatch");
  require(T >= 1 && T <= 5, "FFT3D invalid T");
  require(c >= 0 && c < T, "FFT3D invalid c");
  if (!bins)
    return;

  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  const auto zero = hn::Zero(d), eps = hn::Set(d, 1e-15f);
  const auto noise_v = hn::Set(d, noise.uniform), lower_v = hn::Set(d, lower);
  auto non_finite = hn::MaskFalse(d);

  std::size_t k = 0;
  if (T == 1) {
    const float mr_scale = (grid && degrid != 0.0f && grid[0].real() != 0.0f)
                               ? finite((degrid * spectra[0][0].real()) / grid[0].real())
                               : 0.0f;
    const auto mr_scale_v = hn::Set(d, mr_scale);
    const auto* in_ptr = reinterpret_cast<const float*>(spectra[0]);
    const auto* grid_ptr = reinterpret_cast<const float*>(grid);
    auto* out_ptr = reinterpret_cast<float*>(out);
    for (; bins - k >= lanes; k += lanes) {
      auto re = zero, im = zero, mr = zero, mi = zero;
      hn::LoadInterleaved2(d, in_ptr + 2 * k, re, im);
      if (grid && mr_scale != 0.0f) {
        hn::LoadInterleaved2(d, grid_ptr + 2 * k, mr, mi);
        mr = hn::Mul(mr, mr_scale_v);
        mi = hn::Mul(mi, mr_scale_v);
        re = hn::Sub(re, mr);
        im = hn::Sub(im, mi);
      }
      const auto power = hn::Add(hn::Mul(re, re), hn::Mul(im, im));
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(power)));
      const auto q = hn::Add(power, eps);
      const auto gain = hn::Max(hn::Div(hn::Sub(q, noise.mode == PrimaryMode::Table ? hn::Mul(hn::LoadU(d, noise.table.data() + k), hn::Set(d, noise.multiplier)) : noise_v), q), lower_v);
      non_finite = hn::Or(non_finite, hn::Not(hn::IsFinite(gain)));
      if (grid && mr_scale != 0.0f) {
        re = hn::Add(hn::Mul(gain, re), mr);
        im = hn::Add(hn::Mul(gain, im), mi);
      } else {
        re = hn::Mul(gain, re);
        im = hn::Mul(gain, im);
      }
      non_finite = hn::Or(non_finite, hn::Not(hn::And(hn::IsFinite(re), hn::IsFinite(im))));
      hn::StoreInterleaved2(re, im, d, out_ptr + 2 * k);
    }
  } else {
    const auto& twiddles = get_fft3d_twiddles();
    const float ratio = (grid && degrid != 0 && grid[0].real() != 0)
                          ? finite(finite(degrid*spectra[c][0].real())/grid[0].real()) : 0;
    // Accumulate each temporal frequency directly into the selective inverse.
    // No arrays of sizeless vectors: this also compiles for SVE/RVV targets.
    for (; bins-k >= lanes; k+=lanes) {
      auto grid_re=zero,grid_im=zero;
      if (grid && ratio != 0) {
        hn::LoadInterleaved2(d,reinterpret_cast<const float*>(grid)+2*k,grid_re,grid_im);
        grid_re=hn::Mul(hn::Mul(grid_re,hn::Set(d,ratio)),hn::Set(d,float(T)));
        grid_im=hn::Mul(hn::Mul(grid_im,hn::Set(d,ratio)),hn::Set(d,float(T)));
      }
      auto y_re=zero,y_im=zero;
      for(int m=0;m<T;++m) {
        auto re=zero,im=zero;
        for(int j=0;j<T;++j) {
          auto in_re=zero,in_im=zero;
          hn::LoadInterleaved2(d,reinterpret_cast<const float*>(spectra[j])+2*k,in_re,in_im);
          const auto wr=hn::Set(d,twiddles.fwd[T][m][j].real()),wi=hn::Set(d,twiddles.fwd[T][m][j].imag());
          re=hn::Add(re,hn::Sub(hn::Mul(in_re,wr),hn::Mul(in_im,wi)));
          im=hn::Add(im,hn::Add(hn::Mul(in_im,wr),hn::Mul(in_re,wi)));
        }
        if(m==0) {re=hn::Sub(re,grid_re);im=hn::Sub(im,grid_im);}
        const auto power=hn::Add(hn::Mul(re,re),hn::Mul(im,im));
        non_finite=hn::Or(non_finite,hn::Not(hn::IsFinite(power)));
        const auto q=hn::Add(power,eps);
        const auto power_noise=noise.mode==PrimaryMode::Table
          ? hn::Mul(hn::LoadU(d,noise.table.data()+k),hn::Set(d,noise.multiplier)) : noise_v;
        const auto gain=hn::Max(hn::Div(hn::Sub(q,power_noise),q),lower_v);
        re=hn::Mul(re,gain);im=hn::Mul(im,gain);
        if(m==0) {re=hn::Add(re,grid_re);im=hn::Add(im,grid_im);}
        const auto wr=hn::Set(d,twiddles.inv[T][c][m].real()),wi=hn::Set(d,twiddles.inv[T][c][m].imag());
        y_re=hn::Add(y_re,hn::Sub(hn::Mul(re,wr),hn::Mul(im,wi)));
        y_im=hn::Add(y_im,hn::Add(hn::Mul(im,wr),hn::Mul(re,wi)));
      }
      y_re=hn::Mul(y_re,hn::Set(d,1.0f/float(T)));y_im=hn::Mul(y_im,hn::Set(d,1.0f/float(T)));
      non_finite=hn::Or(non_finite,hn::Not(hn::And(hn::IsFinite(y_re),hn::IsFinite(y_im))));
      hn::StoreInterleaved2(y_re,y_im,d,reinterpret_cast<float*>(out)+2*k);
    }
  }

  if (!hn::AllFalse(d, non_finite))
    throw std::runtime_error("non-finite sample or intermediate");

  if (k < bins) {
    if (T == 1) {
      const float mr_scale = (grid && degrid != 0.0f && grid[0].real() != 0.0f)
                                 ? finite((degrid * spectra[0][0].real()) / grid[0].real())
                                 : 0.0f;
      for (; k < bins; ++k) {
        const float mr = grid ? finite(mr_scale * grid[k].real()) : 0.0f;
        const float mi = grid ? finite(mr_scale * grid[k].imag()) : 0.0f;
        const float re = finite(spectra[0][k].real() - mr);
        const float im = finite(spectra[0][k].imag() - mi);
        const float power = finite(re * re + im * im);
        const float q = power + 1e-15f;
        const float gain = std::max((q - noise.at(k)) / q, lower);
        out[k] = {finite(gain * re + mr), finite(gain * im + mi)};
      }
    } else {
      const float g_ratio = (grid && degrid != 0.0f && grid[0].real() != 0.0f)
                                ? finite((degrid * spectra[c][0].real()) / grid[0].real())
                                : 0.0f;
      const float inv_T = 1.0f / float(T);
      const auto& twiddles = get_fft3d_twiddles();
      const auto* fwd_twiddle = twiddles.fwd[T];
      const auto* inv_twiddle = twiddles.inv[T][c];
      for (; k < bins; ++k) {
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
  }
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
HWY_EXPORT(Fft3dTemporal);
HWY_EXPORT(Target);
HWY_EXPORT(SimdLanes);
SpectralKernel select_spectral(int opt) {
  return opt == 1 ? spectral_scalar : HWY_DYNAMIC_DISPATCH(Spectral);
}
const char* spectral_target(int opt) {
  select_spectral(opt);
  return opt == 1 ? "scalar" : HWY_DYNAMIC_DISPATCH(Target)();
}
Fft3dTemporalKernel select_fft3d_temporal(int opt) {
  return opt == 1 ? fft3d_temporal_scalar : HWY_DYNAMIC_DISPATCH(Fft3dTemporal);
}
const char* fft3d_temporal_target(int opt) {
  select_fft3d_temporal(opt);
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
