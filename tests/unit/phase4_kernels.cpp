#include "kernels/kalman.hpp"
#include "kernels/rows.hpp"
#include "kernels/dither.hpp"
#include "kernels/spatial.hpp"
#include "../guard.hpp"
#include <cstring>
#include <algorithm>
#ifdef NEO_FFT_TEST_HIGHWAY
#include <hwy/targets.h>
#endif
using namespace neo_fft;
using Z=std::complex<float>;
void run_tests() {
  const auto kernel=select_kalman(0);const auto noise=select_dither_noise(0);const auto copy=select_copy_row(0);
  kernel(nullptr,nullptr,nullptr,nullptr,nullptr,0,0,0);noise(nullptr,0,0,0,0,0,0);copy(nullptr,nullptr,0);
  for(std::size_t n=1;n<=257;++n) {
    Guard<std::uint8_t> src(n),dst(n);for(std::size_t i=0;i<n;++i)src.data[i]=std::uint8_t(i*37);
    copy(src.data,dst.data,n);CHECK(std::memcmp(src.data,dst.data,n)==0);
    Guard<float> random(n),expected(n),pattern(n);
    for(auto seed:{0u,17u,2147483647u}) {
      noise(random.data,n,seed,2147483646u,2,107,0);
      dither_noise_scalar(expected.data,n,seed,2147483646u,2,107,0);
      CHECK(std::memcmp(random.data,expected.data,n*sizeof(float))==0);
    }
    Guard<Z> x(n),a(n),b(n),ac(n),bc(n),aq(n),bq(n);
    for(bool table:{false,true}) for(float uniform:{0.f,1.f,2.f}) for(float ratio:{0.f,4.f,100.f}) {
      for(std::size_t i=0;i<n;++i) {
        x.data[i]={float(int(i%7)-3)*.37f,float(int(i%13)-6)*.17f};
        a.data[i]=b.data[i]={float(i%3)*.2f,float(i%5)*.1f};
        ac.data[i]=bc.data[i]={2.f+float(i%4),3.f};aq.data[i]=bq.data[i]={1.f,4.f};
        pattern.data[i]=float(i%5)*.25f;
      }
      for(int iteration=0;iteration<3;++iteration) {
        kernel(x.data,a.data,ac.data,aq.data,table ? pattern.data:nullptr,uniform,ratio,n);
        kalman_scalar(x.data,b.data,bc.data,bq.data,table ? pattern.data:nullptr,uniform,ratio,n);
        CHECK(std::memcmp(a.data,b.data,n*sizeof(Z))==0);
        CHECK(std::memcmp(ac.data,bc.data,n*sizeof(Z))==0 && std::memcmp(aq.data,bq.data,n*sizeof(Z))==0);
      }
    }
  }
  // Long histories amplify a single wrong motion decision; compare every state
  // component after every step, including values adjacent to the strict threshold.
  for(int signal=0;signal<4;++signal) {
    constexpr int n=33;
    Guard<Z> x(n),l(n),c(n),q(n),sl(n),sc(n),sq(n);
    std::fill_n(c.data,n,Z(1,1));std::fill_n(q.data,n,Z(1,1));
    std::fill_n(sc.data,n,Z(1,1));std::fill_n(sq.data,n,Z(1,1));
    for(int frame=0;frame<512;++frame) {
      for(int i=0;i<n;++i) {
        float v=.25f;
        if(signal==1)v=float(frame)*.002f+float(i)*.001f;
        if(signal==2) {
          const float threshold=sl.data[i].real()+2.f;
          v=i%3==0 ? threshold : std::nextafter(threshold,i%3==1 ? INFINITY : -INFINITY);
        }
        if(signal==3)v=(frame/19)%2 ? 50.f : -50.f;
        x.data[i]={v,signal==2 ? 0.f : v*.5f};
      }
      kernel(x.data,l.data,c.data,q.data,nullptr,1,4,n);
      kalman_scalar(x.data,sl.data,sc.data,sq.data,nullptr,1,4,n);
      CHECK(std::memcmp(l.data,sl.data,n*sizeof(Z))==0);
      CHECK(std::memcmp(c.data,sc.data,n*sizeof(Z))==0 && std::memcmp(q.data,sq.data,n*sizeof(Z))==0);
    }
  }
  for(auto kernel:{select_kalman(0),select_kalman(1)}) {
    const int n=33;Z x[n],l[n],c[n],q[n];float pattern[n];
    for(int bad=0;bad<n;++bad) {
      std::fill_n(x,n,Z(1e30f,0));std::fill_n(l,n,Z(-1e30f,0));std::fill_n(c,n,Z(3e38f,3e38f));std::copy_n(c,n,q);
      kernel(x,l,c,q,nullptr,1,4,n); // Inactive covariance overflow is allowed on reset.
      std::fill_n(x,n,Z());std::fill_n(l,n,Z());std::fill_n(c,n,Z(1,1));std::fill_n(q,n,Z(1,1));
      c[bad]=q[bad]=Z(3e38f,3e38f);rejects([&]{kernel(x,l,c,q,nullptr,1,4,n);});
      std::fill_n(c,n,Z(1,1));std::fill_n(q,n,Z(1,1));x[bad]={NAN,0};
      rejects([&]{kernel(x,l,c,q,nullptr,0,4,n);});x[bad]={0,0};std::fill_n(pattern,n,1.f);pattern[bad]=NAN;
      rejects([&]{kernel(x,l,c,q,pattern,1,4,n);});
    }
  }
  float gold;noise(&gold,1,17,23,2,5,7);CHECK(gold==10401337.f*0x1p-24f);
  std::cout<<"Kalman exact state, inactive-branch safety, coordinate hash and copy guard tails passed\n";
}
int main(){try {
#ifdef NEO_FFT_TEST_HIGHWAY
  const auto targets=hwy::SupportedAndGeneratedTargets();CHECK(!targets.empty());
  for(auto target:targets) {hwy::SetSupportedTargetsForTest(target);CHECK(std::string(spatial_target(0))==hwy::TargetName(target));std::cout<<hwy::TargetName(target)<<": ";run_tests();}
  hwy::SetSupportedTargetsForTest(0);
#else
  run_tests();
#endif
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
