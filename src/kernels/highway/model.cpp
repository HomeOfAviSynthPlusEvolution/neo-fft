#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/model.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/model.hpp"
#include "base/checked.hpp"

HWY_BEFORE_NAMESPACE();
namespace neo_fft { namespace HWY_NAMESPACE {
namespace hn=hwy::HWY_NAMESPACE;
template<class T> void Decode(const void* input,float* out,std::size_t n,float base,float scale) {
  const auto* src=static_cast<const T*>(input);
  const hn::ScalableTag<float> d;
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
  if(i<n) model_scalar().decode(src+i,sample_storage<T>,out+i,n-i,base,scale);
}
void DecodeSamples(const void* src,SampleStorage type,float* out,std::size_t n,float base,float scale) {
  switch(type) {
    case SampleStorage::U8: Decode<std::uint8_t>(src,out,n,base,scale);break;
    case SampleStorage::U16: Decode<std::uint16_t>(src,out,n,base,scale);break;
    case SampleStorage::F32: Decode<float>(src,out,n,base,scale);break;
  }
}
void Window(float* values,const float* weights,std::size_t n,float factor) {
  const hn::ScalableTag<float> d;const auto lanes=hn::Lanes(d);auto bad=hn::MaskFalse(d);
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    auto v=hn::LoadU(d,values+i);bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    v=hn::Mul(v,hn::Set(d,factor));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    v=hn::Mul(v,hn::LoadU(d,weights+i));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    hn::StoreU(v,d,values+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if(i<n) model_scalar().window(values+i,weights+i,n-i,factor);
}
void Power(const std::complex<float>* spectrum,const std::complex<float>* grid,float ratio,
           float* out,std::size_t n,bool accumulate) {
  const hn::ScalableTag<float> d;const auto lanes=hn::Lanes(d);auto bad=hn::MaskFalse(d);
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
  if(i<n) model_scalar().power(spectrum+i,grid ? grid+i : nullptr,ratio,out+i,n-i,accumulate);
}
float Score(const float* power,const float* weights,std::size_t n) {
  const hn::ScalableTag<float> d;const auto lanes=hn::Lanes(d);
  HWY_ALIGN float terms[hn::MaxLanes(d)];
  float sum=0;std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    hn::Store(hn::Mul(hn::LoadU(d,power+i),hn::LoadU(d,weights+i)),d,terms);
    for(std::size_t j=0;j<lanes;++j) sum=finite(sum+finite(terms[j]));
  }
  // Do not add a separately reduced tail: preserve the full left-to-right sum.
  for(;i<n;++i) sum=finite(sum+finite(power[i]*weights[i]));
  return sum;
}
void Scale(float* values,const float* weights,std::size_t n,float factor,float maximum) {
  const hn::ScalableTag<float> d;const auto lanes=hn::Lanes(d);auto bad=hn::MaskFalse(d);
  std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    auto v=hn::Mul(hn::Set(d,factor),hn::LoadU(d,values+i));bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    if(weights) v=hn::Mul(v,hn::LoadU(d,weights+i));
    bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));
    bad=hn::Or(bad,hn::Not(hn::IsFinite(hn::Mul(hn::Set(d,maximum),v))));
    hn::StoreU(v,d,values+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if(i<n) model_scalar().scale(values+i,weights ? weights+i : nullptr,n-i,factor,maximum);
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
