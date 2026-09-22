#include "kernels/table.hpp"
#include "base/checked.hpp"
#include <algorithm>
namespace neo_fft { namespace {
void radius(float* out,const float* x,std::size_t n,float other,float dimensions) {
  for(std::size_t i=0;i<n;++i) out[i]=finite(std::sqrt(finite(finite(other+finite(x[i]*x[i]))/dimensions)));
}
void analytic(float* out,const float* freq,std::size_t n,const float* s,float norm) {
  const float a=std::sqrt(.5f)/4,b=std::sqrt(.5f)/2;
  for(std::size_t i=0;i<n;++i) {
    const float f=freq[i];
    const float v=f<a ? finite(s[3]+finite(finite((s[2]-s[3])*f)/a))
                : f<b ? finite(s[2]+finite(finite((s[1]-s[2])*(f-a))/(b-a)))
                      : finite(s[0]+finite(finite((s[1]-s[0])*(1-f))/(1-b)));
    out[i]=finite(finite(v*v)/norm);
  }
}
void curve(float* out,const float* freq,std::size_t n,const float* positions,const float* values,std::size_t knots,float divisor) {
  for(std::size_t i=0;i<n;++i) {
    const auto* hi=std::lower_bound(positions,positions+knots,freq[i]);
    require(hi!=positions+knots,"DFTTest frequency outside curve");const auto k=std::size_t(hi-positions);
    float value;
    if(*hi==freq[i]) value=values[k];
    else {
      require(k>0,"DFTTest frequency outside curve");
      const float t=finite((freq[i]-positions[k-1])/(*hi-positions[k-1]));
      value=finite(finite(values[k-1]*(1-t))+finite(values[k]*t));
    }
    out[i]=finite(value/divisor);
  }
}
void product(float* out,const float* x,std::size_t n,float factor,float divisor) {
  for(std::size_t i=0;i<n;++i) out[i]=finite(finite(factor*x[i])/divisor);
}
void weights(float* out,const float* x,std::size_t n,float other,float cutoff) {
  for(std::size_t i=0;i<n;++i) {
    const float q=finite(other+finite(x[i]*x[i]));out[i]=finite(q/finite(q+cutoff));
  }
}
void enhancement(float* sharpen,float* halo,const float* x2,std::size_t n,float other,float cutoff,float hr) {
  for(std::size_t i=0;i<n;++i) {
    const float d2=finite(other+x2[i]);
    if(sharpen) sharpen[i]=finite(1-std::exp(finite(-d2/cutoff)));
    if(halo) {
      const float first=finite(finite(finite(-.7f*d2)*hr)*hr);
      const float second=finite(finite(-d2*hr)*hr);
      halo[i]=finite(std::exp(first)-std::exp(second));
    }
  }
}
void divide(float* out,std::size_t n,float divisor) {
  for(std::size_t i=0;i<n;++i) out[i]=finite(out[i]/divisor);
}
float window(float* out,const double* x,std::size_t n,double factor,double norm,float energy) {
  for(std::size_t i=0;i<n;++i) {
    const double v=factor*x[i]*norm;
    require(std::isfinite(v) && std::abs(v)<=std::numeric_limits<float>::max(),"DFTTest window overflow");
    const float h=float(v);out[i]=h;energy+=h*h;
  }
  return energy;
}
float maximum(const float* x,std::size_t n) {
  float result=0;for(std::size_t i=0;i<n;++i)result=std::max(result,x[i]);return result;
}
void check_scale(const float* x,std::size_t n,float factor) {
  for(std::size_t i=0;i<n;++i)finite(factor*x[i]);
}
}
const TableKernels& table_scalar() {
  static const TableKernels kernels{radius,analytic,curve,product,weights,enhancement,divide,window,maximum,check_scale};return kernels;
}
#if !NEO_FFT_ENABLE_HIGHWAY
TableKernels select_table(int) {return table_scalar();}
#endif
}
