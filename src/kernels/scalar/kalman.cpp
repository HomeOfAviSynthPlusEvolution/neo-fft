#include "kernels/kalman.hpp"
namespace neo_fft {
void kalman_scalar(const std::complex<float>* x, std::complex<float>* last, std::complex<float>* covariance,
                   std::complex<float>* process, const float* pattern, float uniform, float ratio2, std::size_t count) {
  for (std::size_t k=0;k<count;++k) {
    const float r = pattern ? std::max(finite(pattern[k]),1e-15f) : finite(uniform);
    const float threshold = finite(r*ratio2);
    const float xr=finite(x[k].real()), xi=finite(x[k].imag());
    const float lr=finite(last[k].real()), li=finite(last[k].imag());
    const float cr=finite(covariance[k].real()), ci=finite(covariance[k].imag());
    const float qr=finite(process[k].real()), qi=finite(process[k].imag());
    const float dr=xr-lr, di=xi-li;
    if ((!pattern && r==0) || dr*dr>threshold || di*di>threshold) {
      last[k]=x[k]; covariance[k]=process[k]={r,r};
      continue;
    }
    const auto step = [r](float v,float l,float c,float q,float& nl,float& nc,float& nq) {
      const float sum=finite(c+q);
      const float gain=finite(sum/finite(sum+r));
      nq=finite((gain*gain)*r);
      nc=finite((1-gain)*sum);
      nl=finite(finite(gain*v)+finite((1-gain)*l));
    };
    float nr,ni,ncr,nci,nqr,nqi;
    step(xr,lr,cr,qr,nr,ncr,nqr); step(xi,li,ci,qi,ni,nci,nqi);
    last[k]={nr,ni}; covariance[k]={ncr,nci}; process[k]={nqr,nqi};
  }
}
} // namespace neo_fft
