#include "plugin/filter.hpp"
#include "../test.hpp"
#include <atomic>
#include <thread>
using namespace neo_fft;
int main() {try {
  runtime::Executor executor(3);
  std::atomic<int> entered{0};
  executor.run(3,[&](int) {
    entered.fetch_add(1);
    const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
    while(entered.load()<3 && std::chrono::steady_clock::now()<deadline) std::this_thread::yield();
    CHECK(entered==3);
  });CHECK(executor.peak()==3);
  rejects([&]{executor.run(24,[](int n){if(n==3)throw std::runtime_error("worker failure");});});
  std::atomic<int> done{0};executor.run(24,[&](int){++done;});CHECK(done==24);
  std::vector<float> input(192*128),a(input.size()),b(input.size());
  for(std::size_t i=0;i<input.size();++i)input[i]=float((i*37)%251)/255;
  span2d::Plane<const float> source(input.data(),192,128,192*4);
  span2d::Plane<float> dst_a(a.data(),192,128,192*4),dst_b(b.data(),192,128,192*4);
  for(int T:{1,3}) for(int mode:{0,1}) {
    DFTConfig c;c.block=mode ? 8 : 9;c.overlap=4;c.mode=mode;c.tbsize=T;c.opt=1;c.threads=1;
    Plan serial(192,128,{32,true,false},c);c.threads=3;Plan parallel(192,128,{32,true,false},c);
    std::array<span2d::Plane<const float>,3> sources{source,source,source};
    serial.process({sources.data(),std::size_t(T)},dst_a);parallel.process({sources.data(),std::size_t(T)},dst_b);
    CHECK(a==b);CHECK(parallel.executor().peak()>1 && parallel.executor().peak()<=3);
  }
  {
    using F=plugin::Filter<Algorithm::FFT3D>;
    F::State state;state.source={192,128,1,{ds::ColorFamily::Rgb,ds::SampleFormat::Float32,3,0,0}};
    state.executor=std::make_shared<runtime::Executor>(3);state.retention=std::make_shared<runtime::Retention>(3);
    FFT3DConfig c;c.bw=c.bh=8;c.opt=1;c.mt=true;
    struct Provider : ds::VideoFrameProvider {
      ds::VideoFrameView view;
      ds::Result<ds::RequestedVideoFrame> get(int input,int n) override {return ds::Result<ds::RequestedVideoFrame>::success({input,n,view,{}});}
    } provider;
    provider.view.format=state.source.format;provider.view.plane_count=3;
    ds::MutableVideoFrameView output;output.format=state.source.format;output.plane_count=3;
    std::array<std::vector<float>,3> pixels;
    for(int p=0;p<3;++p) {
      state.plans[p]=std::make_shared<Plan>(192,128,SampleFormat{32,true,false},c,state.executor,state.retention);
      state.rois[p]={0,0,192,128,false};provider.view.planes[p]={input.data(),192*4,192,128};
      pixels[p].resize(input.size());output.planes[p]={pixels[p].data(),192*4,192,128};
    }
    ds::VideoProcessContext ctx{0,provider,output,&state};F::process(ctx);
    CHECK(state.executor->peak()>1 && state.executor->peak()<=3);CHECK(state.retention->count()<=3);
  }
  {DFTConfig c;c.block=4;c.overlap=0;c.tbsize=3;c.swin=6;c.twin=7;c.zmean=true;c.threads=3;
    // Existing admission rejects this degenerate template before processing.
    rejects([&]{Plan p(192,128,{32,true,false},c);});
  }
  const auto geom=Geometry{fft3d_axis(32,8,4),fft3d_axis(24,8,4)};
  const auto budget=runtime::make_workspace_budget(geom,64,40,true,1,1);
  auto retained=std::make_shared<runtime::Retention>(1,64*1024*1024);
  {runtime::WorkspacePool x(budget,1,retained),y(budget,1,retained);
    {auto one=x.acquire(),two=x.acquire(),three=y.acquire();CHECK(x.active_count()==2);}
    CHECK(retained->count()==1);CHECK(x.idle_count()+y.idle_count()==1);
  }CHECK(retained->count()==0 && retained->bytes()==0);
  auto tiny=std::make_shared<runtime::Retention>(1,1);
  {runtime::WorkspacePool p(budget,1,tiny);{auto lease=p.acquire();}CHECK(p.idle_count()==0);}
  std::cout<<"Execution: real parallel work, canonical batch reduction, joined failures and shared idle bounds passed\n";
  return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
