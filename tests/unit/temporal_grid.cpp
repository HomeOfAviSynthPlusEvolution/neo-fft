#include "algorithms/temporal_grid.hpp"
#include "algorithms/windows.hpp"
#include "../test.hpp"
#include <set>
using namespace neo_fft;
int main(){try {
  CHECK(floor_div(-7,3)==-3 && floor_div(-6,3)==-2 && floor_div(-1,3)==-1 && floor_div(7,3)==2);
  for(int T=1;T<=15;++T)for(int O=0;O<T;++O) {
    if(O>T/2 && T%(T-O)) {rejects([&]{temporal_blocks(0,T,T,O);});continue;}
    for(int N:{T,T+1,3*T+7})for(int n=0;n<N;++n) {
      const auto blocks=temporal_blocks(n,N,T,O);std::vector<int> expected;
      // Independent enumeration of the absolute lattice, for small clips only.
      for(int s=-std::max(T-O,O);s<=n;s+=T-O)if(n<s+T)expected.push_back(s);
      CHECK(blocks.size()==expected.size());std::set<int> sources;
      for(std::size_t i=0;i<blocks.size();++i) {
        const auto& b=blocks[i];CHECK(b.start==expected[i] && b.target==n-b.start);
        for(int z=0;z<T;++z) {CHECK(b.slots[z]==std::clamp(expected[i]+z,0,N-1));sources.insert(b.slots[z]);}
        CHECK(b.slots[b.target]==n);
      }
      CHECK(sources.size()<=std::size_t(std::min(N,2*T-1)));
    }
    const auto far=temporal_blocks(INT32_MAX-1,INT32_MAX,T,O);
    CHECK(!far.empty() && far.size()<=std::size_t((T+(T-O)-1)/(T-O)));
    for(const auto& b:far)CHECK(b.target>=0 && b.target<T && b.slots[b.target]==INT32_MAX-1);
    for(int window=0;window<12;++window) {
      const auto w=dft_window_3d(T,1,0,0,7,window,2.5f,2.5f,1,1,O);
      std::vector<double> raw(T);for(int z=0;z<T;++z)raw[z]=dft_raw_window(window,z,T,2.5f);
      float energy=0;
      for(int z=0;z<T;++z) {
        double denominator=0;
        for(int h=z;h>=0;h-=T-O)denominator+=raw[h]*raw[h];
        for(int h=z+T-O;h<T;h+=T-O)denominator+=raw[h]*raw[h];
        const float expected=float((raw[z]/std::sqrt(denominator))*(1/std::sqrt(double(T))));
        CHECK(w.h[z]==expected);energy+=expected*expected;
      }
      CHECK(w.wscale==1.f/energy);
      for(int n=0;n<2*T;++n) {
        float weight=0;for(const auto& b:temporal_blocks(n,2*T,T,O))weight+=float(T)*w.h[b.target]*w.h[b.target];
        CHECK(std::abs(weight-1.f)<5e-7f);
      }
    }
  }
  rejects([]{validate_temporal(4,0,0);});validate_temporal(3,0,INT32_MIN);
  rejects([]{temporal_blocks(0,4,4,4);});rejects([]{temporal_blocks(-1,4,4,0);});
  std::cout<<"Temporal lattice, clamps, bounded far index and all normalized windows passed\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
