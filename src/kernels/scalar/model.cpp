#include "kernels/model.hpp"
#include "base/checked.hpp"

namespace neo_fft {
namespace {
template<class T> void decode(const void* input,float* out,std::size_t n,float base,float scale) {
  const auto* src=static_cast<const T*>(input);
  for(std::size_t i=0;i<n;++i) out[i]=finite(finite(finite(float(src[i]))-base)*scale);
}
void decode_samples(const void* src,SampleStorage storage,float* out,std::size_t n,float base,float scale) {
  switch(storage) {
    case SampleStorage::U8: decode<std::uint8_t>(src,out,n,base,scale); break;
    case SampleStorage::U16: decode<std::uint16_t>(src,out,n,base,scale); break;
    case SampleStorage::F32: decode<float>(src,out,n,base,scale); break;
  }
}
void window(float* values,const float* weights,std::size_t n,float factor) {
  for(std::size_t i=0;i<n;++i) values[i]=finite(finite(finite(values[i])*factor)*weights[i]);
}
void power(const std::complex<float>* spectrum,const std::complex<float>* grid,float ratio,
           float* out,std::size_t n,bool accumulate) {
  for(std::size_t i=0;i<n;++i) {
    const float re=finite(spectrum[i].real()-(grid ? finite(ratio*grid[i].real()) : 0));
    const float im=finite(spectrum[i].imag()-(grid ? finite(ratio*grid[i].imag()) : 0));
    const float q=finite(finite(re*re)+finite(im*im));
    out[i]=accumulate ? finite(out[i]+q) : q;
  }
}
float score(const float* power,const float* weights,std::size_t n) {
  float sum=0;
  for(std::size_t i=0;i<n;++i) sum=finite(sum+finite(power[i]*weights[i]));
  return sum;
}
void scale(float* values,const float* weights,std::size_t n,float factor,float maximum) {
  for(std::size_t i=0;i<n;++i) {
    float v=finite(factor*values[i]);
    if(weights) v=finite(v*weights[i]);
    finite(maximum*v);
    values[i]=v;
  }
}
}
const ModelKernels& model_scalar() {
  static const ModelKernels kernels{decode_samples,window,power,score,scale};
  return kernels;
}
#if !NEO_FFT_ENABLE_HIGHWAY
ModelKernels select_model(int) { return model_scalar(); }
#endif
} // namespace neo_fft
