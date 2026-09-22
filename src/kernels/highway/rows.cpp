#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/rows.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/rows.hpp"
HWY_BEFORE_NAMESPACE();
namespace neo_fft {namespace HWY_NAMESPACE {
namespace hn=hwy::HWY_NAMESPACE;
template<class D> void CopyLanes(const std::uint8_t* src,std::uint8_t* dst,std::size_t n,D d) {
  const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes)hn::StoreU(hn::LoadU(d,src+i),d,dst+i);
  if constexpr(hn::MaxLanes(D{})>1) if(i<n)CopyLanes(src+i,dst+i,n-i,hn::Half<D>{});
}
void Copy(const void* src,void* dst,std::size_t n) {
  CopyLanes(static_cast<const std::uint8_t*>(src),static_cast<std::uint8_t*>(dst),n,hn::ScalableTag<std::uint8_t>{});
}
template<class D> auto Mix(hn::Vec<D> v,D d) {
  v=hn::Xor(v,hn::ShiftRight<16>(v));v=hn::Mul(v,hn::Set(d,0x7feb352du));
  v=hn::Xor(v,hn::ShiftRight<15>(v));v=hn::Mul(v,hn::Set(d,0x846ca68bu));
  return hn::Xor(v,hn::ShiftRight<16>(v));
}
template<class D> void NoiseLanes(float* dst,std::size_t n,std::uint32_t seed,std::uint32_t frame,std::uint32_t plane,std::uint32_t y,std::uint32_t x,D d) {
  const hn::Rebind<std::uint32_t,D> du;const auto lanes=hn::Lanes(d);std::size_t i=0;
  auto h=Mix(hn::Set(du,seed^0xa511e9b3u),du);
  for(const auto v:{frame,plane,y})h=Mix(hn::Xor(h,hn::Set(du,v)),du);
  for(;n-i>=lanes;i+=lanes) {
    const auto hash=Mix(hn::Xor(h,hn::Iota(du,x+std::uint32_t(i))),du);
    const auto u=hn::Mul(hn::ConvertTo(d,hn::ShiftRight<8>(hash)),hn::Set(d,0x1p-24f));
    hn::StoreU(u,d,dst+i);
  }
  if constexpr(hn::MaxLanes(D{})>1) if(i<n)NoiseLanes(dst+i,n-i,seed,frame,plane,y,x+std::uint32_t(i),hn::Half<D>{});
}
void Noise(float* dst,std::size_t n,std::uint32_t seed,std::uint32_t frame,std::uint32_t plane,std::uint32_t y,std::uint32_t x) {
  NoiseLanes(dst,n,seed,frame,plane,y,x,hn::ScalableTag<float>{});
}
}}
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(Copy);HWY_EXPORT(Noise);
CopyRow select_copy_row(int opt) {return opt==1 ? copy_row_scalar : HWY_DYNAMIC_DISPATCH(Copy);}
DitherNoise select_dither_noise(int opt) {return opt==1 ? dither_noise_scalar : HWY_DYNAMIC_DISPATCH(Noise);}
}
#endif
