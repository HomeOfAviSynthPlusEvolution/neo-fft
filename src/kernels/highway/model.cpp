#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/model.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/model.hpp"
#include "base/checked.hpp"

HWY_BEFORE_NAMESPACE();
namespace neo_fft { namespace HWY_NAMESPACE {
namespace hn=hwy::HWY_NAMESPACE;
// Process full vectors first, then narrower exact-width vectors. This keeps
// short sample rows and tails in Highway without masked/padded memory access.
template<class T,class D> void Decode(const void* input,float* out,std::size_t n,float base,float scale,D d) {
  const auto* src=static_cast<const T*>(input);
  const auto lanes=hn::Lanes(d);
  auto bad=hn::MaskFalse(d);
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    auto v=hn::Zero(d);
    if constexpr(std::is_same_v<T,float>) v=hn::LoadU(d,src+i);
    else {
      const hn::Rebind<T,decltype(d)> dt;
      const hn::Rebind<std::uint32_t,decltype(d)> du;
      v=hn::ConvertTo(d,hn::PromoteTo(du,hn::LoadU(dt,src+i)));
    }
    bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    v=hn::Sub(v,hn::Set(d,base));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    v=hn::Mul(v,hn::Set(d,scale));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    hn::StoreU(v,d,out+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if constexpr(hn::MaxLanes(D{})>1)
    if(i<n) Decode<T>(src+i,out+i,n-i,base,scale,hn::Half<D>{});
}
void DecodeSamples(const void* src,SampleStorage type,float* out,std::size_t n,float base,float scale) {
  switch(type) {
    case SampleStorage::U8: Decode<std::uint8_t>(src,out,n,base,scale,hn::ScalableTag<float>{});break;
    case SampleStorage::U16: Decode<std::uint16_t>(src,out,n,base,scale,hn::ScalableTag<float>{});break;
    case SampleStorage::F32: Decode<float>(src,out,n,base,scale,hn::ScalableTag<float>{});break;
  }
}
template<class D> void WindowLanes(float* values,const float* weights,std::size_t n,float factor,D d) {
  const auto lanes=hn::Lanes(d);auto bad=hn::MaskFalse(d);
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    auto v=hn::LoadU(d,values+i);bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    v=hn::Mul(v,hn::Set(d,factor));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    v=hn::Mul(v,hn::LoadU(d,weights+i));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    hn::StoreU(v,d,values+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if constexpr(hn::MaxLanes(D{})>1)
    if(i<n) WindowLanes(values+i,weights+i,n-i,factor,hn::Half<D>{});
}
template<class D> void PowerLanes(const std::complex<float>* spectrum,const std::complex<float>* grid,float ratio,
           float* out,std::size_t n,bool accumulate,D d) {
  const auto lanes=hn::Lanes(d);auto bad=hn::MaskFalse(d);
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    auto re=hn::Zero(d),im=re,mr=re,mi=re;
    hn::LoadInterleaved2(d,reinterpret_cast<const float*>(spectrum)+2*i,re,im);
    if(grid) {
      hn::LoadInterleaved2(d,reinterpret_cast<const float*>(grid)+2*i,mr,mi);
      mr=hn::Mul(mr,hn::Set(d,ratio));mi=hn::Mul(mi,hn::Set(d,ratio));
      bad=hn::Or(bad,hn::Not(hn::And(hn::IsFinite(mr),hn::IsFinite(mi))));
    }
    re=hn::Sub(re,mr);im=hn::Sub(im,mi);
    bad=hn::Or(bad,hn::Not(hn::And(hn::IsFinite(re),hn::IsFinite(im))));
    re=hn::Mul(re,re);im=hn::Mul(im,im);
    bad=hn::Or(bad,hn::Not(hn::And(hn::IsFinite(re),hn::IsFinite(im))));
    auto q=hn::Add(re,im);bad=hn::Or(bad,hn::Not(hn::IsFinite(q)));
    if(accumulate) q=hn::Add(hn::LoadU(d,out+i),q);
    bad=hn::Or(bad,hn::Not(hn::IsFinite(q)));hn::StoreU(q,d,out+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if constexpr(hn::MaxLanes(D{})>1)
    if(i<n) PowerLanes(spectrum+i,grid ? grid+i : nullptr,ratio,out+i,n-i,accumulate,hn::Half<D>{});
}
template<class D> float ScoreLanes(const float* power,const float* weights,std::size_t n,float sum,D d) {
  const auto lanes=hn::Lanes(d);
  HWY_ALIGN float terms[hn::MaxLanes(d)];
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    hn::Store(hn::Mul(hn::LoadU(d,power+i),hn::LoadU(d,weights+i)),d,terms);
    for(std::size_t j=0;j<lanes;++j) sum=finite(sum+finite(terms[j]));
  }
  // Carry the running sum into the narrower width; never sum a tail separately.
  if constexpr(hn::MaxLanes(D{})>1)
    if(i<n) sum=ScoreLanes(power+i,weights+i,n-i,sum,hn::Half<D>{});
  return sum;
}
template<class D> void ScaleLanes(float* values,const float* weights,std::size_t n,float factor,float maximum,D d) {
  const auto lanes=hn::Lanes(d);auto bad=hn::MaskFalse(d);
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    auto v=hn::Mul(hn::Set(d,factor),hn::LoadU(d,values+i));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    if(weights) v=hn::Mul(v,hn::LoadU(d,weights+i));
    bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    bad=hn::Or(bad,hn::Not(hn::IsFinite(hn::Mul(hn::Set(d,maximum),v))));
    hn::StoreU(v,d,values+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if constexpr(hn::MaxLanes(D{})>1)
    if(i<n) ScaleLanes(values+i,weights ? weights+i : nullptr,n-i,factor,maximum,hn::Half<D>{});
}
void Window(float* values,const float* weights,std::size_t n,float factor) {
  WindowLanes(values,weights,n,factor,hn::ScalableTag<float>{});
}
void Power(const std::complex<float>* spectrum,const std::complex<float>* grid,float ratio,
           float* out,std::size_t n,bool accumulate) {
  PowerLanes(spectrum,grid,ratio,out,n,accumulate,hn::ScalableTag<float>{});
}
float Score(const float* power,const float* weights,std::size_t n) {
  return ScoreLanes(power,weights,n,0,hn::ScalableTag<float>{});
}
void Scale(float* values,const float* weights,std::size_t n,float factor,float maximum) {
  ScaleLanes(values,weights,n,factor,maximum,hn::ScalableTag<float>{});
}
ModelKernels GetModelKernels() {return {DecodeSamples,Window,Power,Score,Scale};}
} } // namespace neo_fft::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(GetModelKernels);
ModelKernels select_model(int opt) { return opt==1 ? model_scalar() : HWY_DYNAMIC_DISPATCH(GetModelKernels)(); }
}
#endif
