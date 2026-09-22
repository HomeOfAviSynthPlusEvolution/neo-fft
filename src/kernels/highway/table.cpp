#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/table.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/table.hpp"
#include "base/checked.hpp"
#include <algorithm>
HWY_BEFORE_NAMESPACE();
namespace neo_fft { namespace HWY_NAMESPACE {
namespace hn=hwy::HWY_NAMESPACE;
// Half-spectrum rows are often only 9 or 17 entries; eight lanes avoid turning
// the common small row into an entirely scalar tail on AVX-512 machines.
using D=hn::CappedTag<float,8>;
template<class V> HWY_INLINE V Checked(V v) {
  const D d;
  if(!hn::AllTrue(d,hn::IsFinite(v))) throw std::runtime_error("non-finite sample or intermediate");
  return v;
}
template<class V,class M> HWY_INLINE V CheckedIf(V v,M active) {
  const D d;
  if(!hn::AllFalse(d,hn::And(active,hn::Not(hn::IsFinite(v))))) throw std::runtime_error("non-finite sample or intermediate");
  return v;
}
void Radius(float* out,const float* x,std::size_t n,float other,float dimensions) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    const auto vx=hn::LoadU(d,x+i);
    const auto q=Checked(hn::Div(Checked(hn::Add(hn::Set(d,other),Checked(hn::Mul(vx,vx)))),hn::Set(d,dimensions)));
    hn::StoreU(Checked(hn::Sqrt(q)),d,out+i);
  }
  if(i<n) table_scalar().radius(out+i,x+i,n-i,other,dimensions);
}
void Analytic(float* out,const float* freq,std::size_t n,const float* s,float norm) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  const float a=std::sqrt(.5f)/4,b=std::sqrt(.5f)/2;
  for(;n-i>=lanes;i+=lanes) {
    const auto f=hn::LoadU(d,freq+i);const auto m0=hn::Lt(f,hn::Set(d,a));
    const auto m1=hn::And(hn::Not(m0),hn::Lt(f,hn::Set(d,b)));const auto m2=hn::Not(hn::Or(m0,m1));
    const auto p0=CheckedIf(hn::Mul(hn::Set(d,s[2]-s[3]),f),m0);
    const auto p1=CheckedIf(hn::Mul(hn::Set(d,s[1]-s[2]),hn::Sub(f,hn::Set(d,a))),m1);
    const auto p2=CheckedIf(hn::Mul(hn::Set(d,s[1]-s[0]),hn::Sub(hn::Set(d,1),f)),m2);
    const auto v0=CheckedIf(hn::Add(hn::Set(d,s[3]),CheckedIf(hn::Div(p0,hn::Set(d,a)),m0)),m0);
    const auto v1=CheckedIf(hn::Add(hn::Set(d,s[2]),CheckedIf(hn::Div(p1,hn::Set(d,b-a)),m1)),m1);
    const auto v2=CheckedIf(hn::Add(hn::Set(d,s[0]),CheckedIf(hn::Div(p2,hn::Set(d,1-b)),m2)),m2);
    const auto v=hn::IfThenElse(m0,v0,hn::IfThenElse(m1,v1,v2));
    hn::StoreU(Checked(hn::Div(Checked(hn::Mul(v,v)),hn::Set(d,norm))),d,out+i);
  }
  if(i<n) table_scalar().analytic(out+i,freq+i,n-i,s,norm);
}
void Curve(float* out,const float* freq,std::size_t n,const float* positions,const float* values,std::size_t knots,float divisor) {
  // Gather indices are int32. Very large knot arrays keep the scalar size_t path.
  if(knots>std::size_t(INT32_MAX)) {table_scalar().curve(out,freq,n,positions,values,knots,divisor);return;}
  const D d;const hn::Rebind<std::int32_t,D> di;const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    const auto f=hn::LoadU(d,freq+i);
    require(hn::AllTrue(d,hn::And(hn::Ge(f,hn::Set(d,positions[0])),hn::Le(f,hn::Set(d,positions[knots-1])))),"DFTTest frequency outside curve");
    auto lo=hn::Zero(di),hi=hn::Set(di,int(knots-1));
    while(!hn::AllTrue(di,hn::Eq(lo,hi))) {
      const auto mid=hn::Add(lo,hn::ShiftRight<1>(hn::Sub(hi,lo)));
      const auto less=hn::RebindMask(di,hn::Lt(hn::GatherIndex(d,positions,mid),f));
      lo=hn::IfThenElse(less,hn::Add(mid,hn::Set(di,1)),lo);hi=hn::IfThenElse(less,hi,mid);
    }
    const auto hp=hn::GatherIndex(d,positions,hi),hv=hn::GatherIndex(d,values,hi);
    const auto between=hn::Ne(hp,f);const auto previous=hn::Max(hn::Zero(di),hn::Sub(hi,hn::Set(di,1)));
    const auto lp=hn::GatherIndex(d,positions,previous),lv=hn::GatherIndex(d,values,previous);
    const auto den=hn::IfThenElse(between,hn::Sub(hp,lp),hn::Set(d,1));
    const auto t=CheckedIf(hn::Div(hn::Sub(f,lp),den),between);
    const auto v=CheckedIf(hn::Add(CheckedIf(hn::Mul(lv,hn::Sub(hn::Set(d,1),t)),between),CheckedIf(hn::Mul(hv,t),between)),between);
    hn::StoreU(Checked(hn::Div(hn::IfThenElse(between,v,hv),hn::Set(d,divisor))),d,out+i);
  }
  if(i<n) table_scalar().curve(out+i,freq+i,n-i,positions,values,knots,divisor);
}
void Product(float* out,const float* x,std::size_t n,float factor,float divisor) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes)
    hn::StoreU(Checked(hn::Div(Checked(hn::Mul(hn::Set(d,factor),hn::LoadU(d,x+i))),hn::Set(d,divisor))),d,out+i);
  if(i<n) table_scalar().product(out+i,x+i,n-i,factor,divisor);
}
void Weights(float* out,const float* x,std::size_t n,float other,float cutoff) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) {
    const auto vx=hn::LoadU(d,x+i);const auto q=Checked(hn::Add(hn::Set(d,other),Checked(hn::Mul(vx,vx))));
    hn::StoreU(Checked(hn::Div(q,Checked(hn::Add(q,hn::Set(d,cutoff))))),d,out+i);
  }
  if(i<n) table_scalar().weights(out+i,x+i,n-i,other,cutoff);
}
void Enhancement(float* sharpen,float* halo,const float* x2,std::size_t n,float other,float cutoff,float hr) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  HWY_ALIGN float first[hn::MaxLanes(d)],second[hn::MaxLanes(d)];
  for(;n-i>=lanes;i+=lanes) {
    const auto d2=Checked(hn::Add(hn::Set(d,other),hn::LoadU(d,x2+i)));
    if(sharpen) {
      hn::Store(Checked(hn::Div(hn::Neg(d2),hn::Set(d,cutoff))),d,first);
      for(std::size_t j=0;j<lanes;++j) sharpen[i+j]=finite(1-std::exp(first[j]));
    }
    if(halo) {
      hn::Store(Checked(hn::Mul(Checked(hn::Mul(Checked(hn::Mul(hn::Set(d,-.7f),d2)),hn::Set(d,hr))),hn::Set(d,hr))),d,first);
      hn::Store(Checked(hn::Mul(Checked(hn::Mul(hn::Neg(d2),hn::Set(d,hr))),hn::Set(d,hr))),d,second);
      for(std::size_t j=0;j<lanes;++j) halo[i+j]=finite(std::exp(first[j])-std::exp(second[j]));
    }
  }
  if(i<n) table_scalar().enhancement(sharpen ? sharpen+i : nullptr,halo ? halo+i : nullptr,x2+i,n-i,other,cutoff,hr);
}
void Divide(float* values,std::size_t n,float divisor) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes) hn::StoreU(Checked(hn::Div(hn::LoadU(d,values+i),hn::Set(d,divisor))),d,values+i);
  if(i<n) table_scalar().divide(values+i,n-i,divisor);
}
float Window(float* out,const double* x,std::size_t n,double factor,double norm,float energy) {
  std::size_t i=0;
#if HWY_HAVE_FLOAT64
  const hn::CappedTag<double,4> dd;const hn::Rebind<float,decltype(dd)> df;
  const auto lanes=hn::Lanes(dd);HWY_ALIGN float squares[hn::MaxLanes(df)];
  for(;n-i>=lanes;i+=lanes) {
    const auto v=hn::Mul(hn::Mul(hn::Set(dd,factor),hn::LoadU(dd,x+i)),hn::Set(dd,norm));
    require(hn::AllTrue(dd,hn::And(hn::IsFinite(v),hn::Le(hn::Abs(v),hn::Set(dd,double(std::numeric_limits<float>::max()))))),"DFTTest window overflow");
    const auto h=hn::DemoteTo(df,v);hn::StoreU(h,df,out+i);hn::Store(hn::Mul(h,h),df,squares);
    for(std::size_t j=0;j<lanes;++j)energy+=squares[j];
  }
#endif
  if(i<n)energy=table_scalar().window(out+i,x+i,n-i,factor,norm,energy);
  return energy;
}
float Maximum(const float* x,std::size_t n) {
  const D d;const auto lanes=hn::Lanes(d);auto v=hn::Zero(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes)v=hn::Max(v,hn::LoadU(d,x+i));
  float result=hn::ReduceMax(d,v);
  for(;i<n;++i)result=std::max(result,x[i]);
  return result;
}
void CheckScale(const float* x,std::size_t n,float factor) {
  const D d;const auto lanes=hn::Lanes(d);std::size_t i=0;
  for(;n-i>=lanes;i+=lanes)Checked(hn::Mul(hn::Set(d,factor),hn::LoadU(d,x+i)));
  if(i<n)table_scalar().check_scale(x+i,n-i,factor);
}
TableKernels GetTableKernels() {return {Radius,Analytic,Curve,Product,Weights,Enhancement,Divide,Window,Maximum,CheckScale};}
} }
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(GetTableKernels);
TableKernels select_table(int opt) {return opt==1 ? table_scalar() : HWY_DYNAMIC_DISPATCH(GetTableKernels)();}
}
#endif
