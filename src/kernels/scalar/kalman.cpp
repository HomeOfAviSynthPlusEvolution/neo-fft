#include "kernels/kalman.hpp"
namespace neo_fft {
void kalman_scalar(const std::complex<float>* x, std::complex<float>* last, float* covariance,
                   float* process, const float* pattern, float uniform, float ratio2, std::size_t count) {
  for (std::size_t k=0;k<count;++k) {
    const float r = pattern ? std::max(finite(pattern[k]),1e-15f) : finite(uniform);
    const float threshold = finite(r*ratio2);
    const float xr=finite(x[k].real()), xi=finite(x[k].imag());
    const float lr=finite(last[k].real()), li=finite(last[k].imag());
    const float c=finite(covariance[k]), q=finite(process[k]);
    const float dr=xr-lr, di=xi-li;
    if ((!pattern && r==0) || dr*dr>threshold || di*di>threshold) {
      last[k]=x[k]; covariance[k]=process[k]=r;
      continue;
    }
    const float sum=finite(c+q);
    const float gain=finite(sum/finite(sum+r));
    const float complement=1-gain;
    const float nq=finite((gain*gain)*r), nc=finite(complement*sum);
    const float nr=finite(finite(gain*xr)+finite(complement*lr));
    const float ni=finite(finite(gain*xi)+finite(complement*li));
    last[k]={nr,ni}; covariance[k]=nc; process[k]=nq;
  }
}
#if !NEO_FFT_ENABLE_HIGHWAY
KalmanKernel select_kalman(int) {return kalman_scalar;}
#endif
} // namespace neo_fft
