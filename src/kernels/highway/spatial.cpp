#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/spatial.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/spatial.hpp"
#include "base/checked.hpp"
#include <algorithm>

HWY_BEFORE_NAMESPACE();
namespace neo_fft {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

void ValidateFinite(const float* src, std::size_t count) {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  auto invalid = hn::MaskFalse(d);
  std::size_t i = 0;
  for (; count - i >= lanes; i += lanes)
    invalid = hn::Or(invalid, hn::Not(hn::IsFinite(hn::LoadU(d, src + i))));
  if (!hn::AllFalse(d, invalid))
    throw std::runtime_error("non-finite sample or intermediate");
  for (; i < count; ++i) finite(src[i]);
}

void GatherFft3d(const float* src, const float* wx_a, float wy, float* blk, int count) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  const auto vwy = hn::Set(d, wy);
  int x = 0;
  for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
    const auto vsrc = hn::LoadU(d, src + x);
    const auto vwx = hn::LoadU(d, wx_a + x);
    hn::StoreU(hn::Mul(hn::Mul(vsrc, vwy), vwx), d, blk + x);
  }
  for (; x < count; ++x) {
    blk[x] = (src[x] * wy) * wx_a[x];
  }
}

void GatherDfttest(const float* src, const float* h_row, float* blk, int count) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  int x = 0;
  for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
    const auto vsrc = hn::LoadU(d, src + x);
    const auto vh = hn::LoadU(d, h_row + x);
    hn::StoreU(hn::Mul(vsrc, vh), d, blk + x);
  }
  for (; x < count; ++x) {
    blk[x] = src[x] * h_row[x];
  }
}

void ScatterFft3dBlock(const float* inv, const float* wx_s, float* r_row, int count) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  int x = 0;
  for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
    const auto vinv = hn::LoadU(d, inv + x);
    const auto vwx = hn::LoadU(d, wx_s + x);
    const auto vr = hn::LoadU(d, r_row + x);
    hn::StoreU(hn::MulAdd(vinv, vwx, vr), d, r_row + x);
  }
  for (; x < count; ++x) {
    r_row[x] += inv[x] * wx_s[x];
  }
}

void ScatterFft3dRow(const float* r_ptr, float wy, float* a_ptr, int count) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  const auto vwy = hn::Set(d, wy);
  int x = 0;
  for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
    const auto vr = hn::LoadU(d, r_ptr + x);
    const auto va = hn::LoadU(d, a_ptr + x);
    hn::StoreU(hn::MulAdd(vr, vwy, va), d, a_ptr + x);
  }
  for (; x < count; ++x) {
    a_ptr[x] += r_ptr[x] * wy;
  }
}

void ScatterDfttest(const float* inv, const float* h_syn, float* acc_row, int count) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  int x = 0;
  for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
    const auto vinv = hn::LoadU(d, inv + x);
    const auto vh = hn::LoadU(d, h_syn + x);
    const auto vacc = hn::LoadU(d, acc_row + x);
    hn::StoreU(hn::MulAdd(vinv, vh, vacc), d, acc_row + x);
  }
  for (; x < count; ++x) {
    acc_row[x] += inv[x] * h_syn[x];
  }
}

void StoreOutputFloat(const float* a_ptr, float* dst_row, int count, bool fft3d, float scale) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  int x = 0;
  if (fft3d) {
    const auto zero = hn::Zero(d);
    const auto one = hn::Set(d, 1.0f);
    for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
      const auto va = hn::LoadU(d, a_ptr + x);
      hn::StoreU(hn::Clamp(va, zero, one), d, dst_row + x);
    }
    for (; x < count; ++x) {
      dst_row[x] = std::clamp(a_ptr[x], 0.0f, 1.0f);
    }
  } else {
    const auto vscale = hn::Set(d, scale);
    for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
      const auto va = hn::LoadU(d, a_ptr + x);
      hn::StoreU(hn::Mul(va, vscale), d, dst_row + x);
    }
    for (; x < count; ++x) {
      dst_row[x] = a_ptr[x] * scale;
    }
  }
}

void StoreOutputUint8(const float* a_ptr, std::uint8_t* dst_row, int count, bool fft3d, float base,
                      float scale, float peak) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  const auto zero = hn::Zero(d);
  const auto vpeak = hn::Set(d, peak);
  const auto vhalf = hn::Set(d, 0.5f);
  constexpr hn::ScalableTag<float> df;
  HWY_ALIGN float buf[hn::MaxLanes(df)];
  int x = 0;
  if (fft3d) {
    const auto vbase = hn::Set(d, base);
    for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
      const auto va = hn::LoadU(d, a_ptr + x);
      const auto v = hn::Add(hn::Add(va, vhalf), vbase);
      const auto vc = hn::Clamp(v, zero, vpeak);
      hn::Store(vc, d, buf);
      for (std::size_t i = 0; i < lanes; ++i) {
        dst_row[x + i] = static_cast<std::uint8_t>(buf[i]);
      }
    }
  } else {
    const auto vscale = hn::Set(d, scale);
    for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
      const auto va = hn::LoadU(d, a_ptr + x);
      const auto v = hn::MulAdd(va, vscale, vhalf);
      const auto vc = hn::Clamp(v, zero, vpeak);
      hn::Store(vc, d, buf);
      for (std::size_t i = 0; i < lanes; ++i) {
        dst_row[x + i] = static_cast<std::uint8_t>(buf[i]);
      }
    }
  }
  for (; x < count; ++x) {
    const float v = fft3d ? ((a_ptr[x] + 0.5f) + base) : ((a_ptr[x] * scale) + 0.5f);
    dst_row[x] = static_cast<std::uint8_t>(std::clamp(v, 0.0f, peak));
  }
}

void StoreOutputUint16(const float* a_ptr, std::uint16_t* dst_row, int count, bool fft3d, float base,
                       float scale, float peak) noexcept {
  const hn::ScalableTag<float> d;
  const auto lanes = hn::Lanes(d);
  const auto zero = hn::Zero(d);
  const auto vpeak = hn::Set(d, peak);
  const auto vhalf = hn::Set(d, 0.5f);
  constexpr hn::ScalableTag<float> df;
  HWY_ALIGN float buf[hn::MaxLanes(df)];
  int x = 0;
  if (fft3d) {
    const auto vbase = hn::Set(d, base);
    for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
      const auto va = hn::LoadU(d, a_ptr + x);
      const auto v = hn::Add(hn::Add(va, vhalf), vbase);
      const auto vc = hn::Clamp(v, zero, vpeak);
      hn::Store(vc, d, buf);
      for (std::size_t i = 0; i < lanes; ++i) {
        dst_row[x + i] = static_cast<std::uint16_t>(buf[i]);
      }
    }
  } else {
    const auto vscale = hn::Set(d, scale);
    for (; x + static_cast<int>(lanes) <= count; x += static_cast<int>(lanes)) {
      const auto va = hn::LoadU(d, a_ptr + x);
      const auto v = hn::MulAdd(va, vscale, vhalf);
      const auto vc = hn::Clamp(v, zero, vpeak);
      hn::Store(vc, d, buf);
      for (std::size_t i = 0; i < lanes; ++i) {
        dst_row[x + i] = static_cast<std::uint16_t>(buf[i]);
      }
    }
  }
  for (; x < count; ++x) {
    const float v = fft3d ? ((a_ptr[x] + 0.5f) + base) : ((a_ptr[x] * scale) + 0.5f);
    dst_row[x] = static_cast<std::uint16_t>(std::clamp(v, 0.0f, peak));
  }
}

SpatialKernels GetSpatialKernels() {
  return SpatialKernels{
      ValidateFinite,
      GatherFft3d,
      GatherDfttest,
      ScatterFft3dBlock,
      ScatterFft3dRow,
      ScatterDfttest,
      StoreOutputFloat,
      StoreOutputUint8,
      StoreOutputUint16,
  };
}

const char* SpatialTarget() {
  return hwy::TargetName(HWY_TARGET);
}

} // namespace HWY_NAMESPACE
} // namespace neo_fft
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(GetSpatialKernels);
HWY_EXPORT(SpatialTarget);

SpatialKernels select_spatial(int opt) {
  if (opt == 1) {
    return spatial_scalar();
  }
  return HWY_DYNAMIC_DISPATCH(GetSpatialKernels)();
}

const char* spatial_target(int opt) {
  if (opt == 1) {
    return "scalar";
  }
  return HWY_DYNAMIC_DISPATCH(SpatialTarget)();
}

} // namespace neo_fft
#endif
