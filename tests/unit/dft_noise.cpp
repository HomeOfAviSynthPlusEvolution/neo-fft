#include "algorithms/plan.hpp"
#include "plugin/filter.hpp"
#include "../reference/dft.hpp"
#include "../test.hpp"
#include <future>
#include <set>
using namespace neo_fft;

int main() { try {
  // Independent Hann/rectangular windows and direct DFT, including heavy overlap
  // where h2 must not inherit output normalization. Table tolerance: 0.2% + 0.002.
  for (int T : {1,3,5,15}) for (int S : {1,3,4}) for (bool mean : {false,true})
  for (int window : {0,7}) for (int overlap : {0,S-1}) for (float alpha : {.5f,2.f}) {
    DFTConfig c; c.tbsize=T;c.block=S;c.overlap=overlap;c.swin=window;c.twin=window;c.zmean=mean;c.alpha=alpha;
    c.locations={{0,0,0,0},{2,1,1,2},{0,0,0,0}};
    DFTNoise model(c);
    auto fill=[&](const NoiseLocation& n,int z,span2d::Span<float> out) {
      for(int y=0;y<S;++y) for(int x=0;x<S;++x) out[y*S+x]=float((n.frame+z)*7+n.plane*13+n.y+y+(n.x+x)*3);
    };
    auto first=std::async(std::launch::async,[&]{model.prepare(fill);});
    model.prepare(fill);first.get();auto table=model.power();CHECK(table);
    std::vector<float> h(T*S*S),g(T*S*S),input(T*S*S);
    auto raw=[&](int i,int length) {return window==7 ? 1. : .5-.5*std::cos(2*std::acos(-1.)*(i+.5)/length);};
    double E1=0,E2=0;
    for(int z=0;z<T;++z) for(int y=0;y<S;++y) for(int x=0;x<S;++x) {
      auto axis=[&](int j){double e=0;for(int i=j%(S-overlap);i<S;i+=S-overlap)e+=std::pow(raw(i,S),2);return raw(j,S)/std::sqrt(e);};
      const auto k=(z*S+y)*S+x;
      h[k]=float(raw(z,T)*raw(y,S)*raw(x,S)/std::sqrt(double(T*S*S)));g[k]=255*h[k];
      E2+=h[k]*h[k];E1+=std::pow(raw(z,T)*axis(y)*axis(x)/std::sqrt(double(T*S*S)),2);
    }
    const auto G=direct_dft_3d(g.data(),T,S,S);
    std::vector<double> expected(table->size());
    for(const auto& n:c.locations) {
      for(int z=0;z<T;++z) fill(n,z,{input.data()+z*S*S,std::size_t(S*S)});
      for(std::size_t k=0;k<input.size();++k)input[k]*=h[k];
      const auto X=direct_dft_3d(input.data(),T,S,S);const double ratio=mean ? X[0].real()/G[0].real() : 0;
      for(std::size_t k=0;k<X.size();++k)expected[k]+=std::norm(X[k]-ratio*G[k])*E1/E2*alpha/c.locations.size();
    }
    for(std::size_t k=0;k<expected.size();++k) CHECK(std::abs((*table)[k]-expected[k])<=.002+.002*std::abs(expected[k]));
    model.prepare([](const NoiseLocation&,int,span2d::Span<float>){throw std::runtime_error("warm model must not gather");});
    CHECK(model.power()==table);
  }
  DFTConfig c;c.block=3;c.overlap=1;c.tbsize=3;c.locations={{0,0,0,0}};
  DFTNoise retry(c);
  rejects([&]{retry.prepare([](const NoiseLocation&,int,span2d::Span<float> out){std::fill(out.begin(),out.end(),NAN);});});
  CHECK(!retry.power());
  rejects([&]{retry.prepare([](const NoiseLocation&,int z,span2d::Span<float> out){
    if(z==1)throw std::runtime_error("simulated upstream sample error");
    std::fill(out.begin(),out.end(),10);
  });});
  CHECK(!retry.power());
  retry.prepare([](const NoiseLocation&,int,span2d::Span<float> out){std::fill(out.begin(),out.end(),1);});CHECK(retry.power());
  // Exact declared closure remains identical before and after publication.
  using F=plugin::Filter<Algorithm::DFTTest>;
  F::State state{};state.source.num_frames=20;state.temporal_size=3;
  c.locations={{0,0,0,0},{10,0,0,0},{10,0,0,0}};
  auto model=std::make_shared<DFTNoise>(c);state.dft_noise=model;
  state.plans[0]=std::make_shared<Plan>(12,12,SampleFormat{8,false,false},c,model);
  for(bool warm:{false,true}) {
    if(warm)model->prepare([](const NoiseLocation&,int,span2d::Span<float> out){std::fill(out.begin(),out.end(),0);});
    std::vector<ds::VideoFrameRequest> requests;ds::VideoRequestContext ctx{6,requests,{},&state};F::request(ctx);
    std::set<int> actual;for(auto r:requests)actual.insert(r.frame_number);
    CHECK(actual==std::set<int>({0,1,2,5,6,7,10,11,12}));CHECK(requests.size()==actual.size());
  }
  using FFT=plugin::Filter<Algorithm::FFT3D>;
  FFT3DConfig fc;fc.bw=fc.bh=4;fc.bt=4;fc.pfactor=1;fc.px=fc.py=1;
  FFT::State fs{};fs.source.num_frames=9;fs.temporal_size=4;fs.sampled=true;fs.pattern_frame=8;
  auto fp=std::make_shared<Plan>(20,20,SampleFormat{8,false,false},fc);fs.plans[0]=fp;
  std::vector<std::uint8_t> pixels(400,64);
  for(bool warm:{false,true}) {
    if(warm)fp->prepare_pattern({pixels.data(),20,20,20});
    for(int n:{0,4,8}) {
      std::vector<ds::VideoFrameRequest> requests;ds::VideoRequestContext ctx{n,requests,{},&fs};FFT::request(ctx);
      std::set<int> actual;for(auto r:requests)actual.insert(r.frame_number);
      const auto want=n==4 ? std::set<int>{2,3,4,5,8} : std::set<int>{n,8};
      CHECK(actual==want);CHECK(requests.size()==actual.size());
    }
  }
  std::cout<<"DFT noise direct DFT, independent windows, publication and dependencies passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;} }
