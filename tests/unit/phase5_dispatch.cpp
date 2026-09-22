#include "algorithms/plan.hpp"
#include "../guard.hpp"
#ifdef NEO_FFT_TEST_HIGHWAY
#include <hwy/targets.h>
#endif
using namespace neo_fft;
void verify() {
  constexpr int W=17,H=13,N=17;
  std::vector<std::unique_ptr<Guard<float>>> data;
  for(int n=0;n<N;++n) {
    data.push_back(std::make_unique<Guard<float>>(W*H));
    for(int i=0;i<W*H;++i)data.back()->data[i]=float((n*n*7+i*13)%127)/128.f-.25f;
  }
  Guard<float> actual(W*H),expected(W*H);
  for(int T:{1,2,3,4,5,6,15})for(int O:{0,T/2,T-1})for(int mode:{0,1})for(int ft=0;ft<5;++ft) {
    DFTConfig c;c.tbsize=T;c.temporal_mode=1;c.temporal_overlap=O;c.mode=mode;c.block=mode ? 4 : 5;c.overlap=2;c.ftype=ft;c.sigma=ft<2 ? 8.f : .7f;c.opt=1;
    Plan scalar(W,H,{32,true,false},c);c.opt=0;Plan simd(W,H,{32,true,false},c);
    for(int n:{0,8,16}) {
      std::vector<span2d::Plane<const float>> slots;std::vector<int> targets;
      for(const auto& b:temporal_blocks(n,N,T,O)) {
        targets.push_back(b.target);
        for(int z=0;z<T;++z)slots.emplace_back(data[b.slots[z]]->data,W,H,W*4);
      }
      scalar.process_at<float>({slots.data(),slots.size()},{expected.data,W,H,W*4},n,0,{targets.data(),targets.size()});
      simd.process_at<float>({slots.data(),slots.size()},{actual.data,W,H,W*4},n,0,{targets.data(),targets.size()});
      for(int i=0;i<W*H;++i)CHECK(std::isfinite(actual.data[i]) && std::abs(actual.data[i]-expected.data[i])<=3e-6f);
    }
  }
  for(int T=1;T<=15;++T)for(int O=0;O<T;++O) {
    if(O>T/2 && T%(T-O))continue;
    for(int w=0;w<12;++w) {
      auto a=dft_window_3d(T,5,2,1,w,w,2.5f,2.5f,0,1,O),b=dft_window_3d(T,5,2,1,w,w,2.5f,2.5f,1,1,O);
      CHECK(a.h==b.h && a.wscale==b.wscale);
    }
  }
}
int main(){try {
#ifdef NEO_FFT_TEST_HIGHWAY
  const auto targets=hwy::SupportedAndGeneratedTargets();CHECK(!targets.empty());
  for(auto target:targets) {
    hwy::SetSupportedTargetsForTest(target);CHECK(std::string(spatial_target(0))==hwy::TargetName(target));
    verify();std::cout<<hwy::TargetName(target)<<": temporal windows/operator/guard planes passed\n";
  }
  hwy::SetSupportedTargetsForTest(0);
#else
  verify();
#endif
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
