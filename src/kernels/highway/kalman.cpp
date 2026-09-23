#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/highway/kalman.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "kernels/kalman.hpp"
HWY_BEFORE_NAMESPACE();
namespace neo_fft {namespace HWY_NAMESPACE {
namespace hn=hwy::HWY_NAMESPACE;
// Clang does not inherit the surrounding target attribute into lambdas.
template<class V,class M> struct Checked {
  M& bad;
  HWY_ATTR V operator()(V v) const {bad=hn::Or(bad,hn::Not(hn::IsFinite(v)));return v;}
};
template<class D> void KalmanLanes(const std::complex<float>* x,std::complex<float>* last,float* covariance,
                                  float* process,const float* pattern,float uniform,float ratio2,std::size_t count,D d) {
  const auto lanes=hn::Lanes(d);const auto zero=hn::Zero(d),one=hn::Set(d,1.f);
  auto bad=hn::MaskFalse(d);
  Checked<decltype(zero),decltype(bad)> check{bad};
  std::size_t i=0;
  for(;count-i>=lanes;i+=lanes) {
    auto xr=zero,xi=zero,lr=zero,li=zero;
    hn::LoadInterleaved2(d,reinterpret_cast<const float*>(x+i),xr,xi);
    hn::LoadInterleaved2(d,reinterpret_cast<const float*>(last+i),lr,li);
    const auto c=check(hn::LoadU(d,covariance+i)),q=check(hn::LoadU(d,process+i));
    check(xr);check(xi);check(lr);check(li);
    const auto r=pattern ? hn::Max(check(hn::LoadU(d,pattern+i)),hn::Set(d,1e-15f)) : check(hn::Set(d,uniform));
    const auto threshold=check(hn::Mul(r,hn::Set(d,ratio2)));
    const auto dr=hn::Sub(xr,lr),di=hn::Sub(xi,li);
    const auto reset=(!pattern && uniform==0) ? hn::Eq(zero,zero)
        : hn::Or(hn::Gt(hn::Mul(dr,dr),threshold),hn::Gt(hn::Mul(di,di),threshold));
    // Inactive smoothing lanes must not divide by zero or overflow C+Q.
    const auto smooth_r=hn::IfThenElse(reset,one,r);
    const auto sum=check(hn::Add(hn::IfThenElse(reset,zero,c),hn::IfThenElse(reset,zero,q)));
    const auto gain=check(hn::Div(sum,check(hn::Add(sum,smooth_r))));
    const auto complement=hn::Sub(one,gain);
    const auto nq=hn::IfThenElse(reset,r,check(hn::Mul(hn::Mul(gain,gain),smooth_r)));
    const auto nc=hn::IfThenElse(reset,r,check(hn::Mul(complement,sum)));
    const auto nr=hn::IfThenElse(reset,xr,check(hn::Add(
        check(hn::Mul(gain,hn::IfThenElse(reset,zero,xr))),
        check(hn::Mul(complement,hn::IfThenElse(reset,zero,lr))))));
    const auto ni=hn::IfThenElse(reset,xi,check(hn::Add(
        check(hn::Mul(gain,hn::IfThenElse(reset,zero,xi))),
        check(hn::Mul(complement,hn::IfThenElse(reset,zero,li))))));
    hn::StoreInterleaved2(nr,ni,d,reinterpret_cast<float*>(last+i));
    hn::StoreU(nc,d,covariance+i);
    hn::StoreU(nq,d,process+i);
  }
  if(!hn::AllFalse(d,bad)) throw std::runtime_error("non-finite sample or intermediate");
  if constexpr(hn::MaxLanes(D{})>1)
    if(i<count) KalmanLanes(x+i,last+i,covariance+i,process+i,pattern ? pattern+i : nullptr,uniform,ratio2,count-i,hn::Half<D>{});
}
void Kalman(const std::complex<float>* x,std::complex<float>* last,float* covariance,
            float* process,const float* pattern,float uniform,float ratio2,std::size_t count) {
  KalmanLanes(x,last,covariance,process,pattern,uniform,ratio2,count,hn::ScalableTag<float>{});
}
}}
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_fft {
HWY_EXPORT(Kalman);
KalmanKernel select_kalman(int opt) {return opt==1 ? kalman_scalar : HWY_DYNAMIC_DISPATCH(Kalman);}
}
#endif
