#include "plugin/filter.hpp"
#include "../test.hpp"
#include <future>
using namespace neo_fft;
using Filter=plugin::Filter<Algorithm::FFT3D>;
struct Pixels : ds::FrameStorage {
  std::vector<float> pixels=std::vector<float>(32*24);
  int* live;
  explicit Pixels(int* counter=nullptr):live(counter) {if(live)++*live;}
  ~Pixels() override {if(live)--*live;}
  ds::VideoFrameView read() const override {
    return {{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0},1,{{{pixels.data(),32*4,32,24}}}};
  }
  ds::MutableVideoFrameView write() override {
    return {{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0},1,{{{pixels.data(),32*4,32,24}}}};
  }
};
struct Factory : ds::FrameFactory {
  bool fail=false;
  ds::WritableFrame allocate(ds::VideoFormat,int,int,const ds::FrameRef&) override {
    if(fail) throw std::bad_alloc();
    return ds::WritableFrame(std::make_unique<Pixels>());
  }
};
Filter::State state(FFT3DConfig config, std::size_t budget=64*1024*1024) {
  config.bt=0;config.bw=8;config.bh=8;config.ow=4;config.oh=4;config.opt=1;
  Filter::State s;
  s.source={32,24,INT32_MAX,{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0}};
  s.executor=std::make_shared<runtime::Executor>(1);
  s.plans[0]=std::make_shared<Plan>(32,24,SampleFormat{32,true,false},config);
  s.rois[0]={0,0,32,24,false};s.kalman=s.plans[0]->kalman();
  s.sampled=s.plans[0]->needs_pattern_frame();s.pattern_frame=11;
  s.checkpoints=std::make_shared<runtime::Checkpoints>(budget);return s;
}
struct Trace {std::vector<int> sources;int peak=0;};
std::vector<float> run(Filter::State& s,int n,Trace& trace,int fail_source=-1,bool fail_output=false) {
  int live=0;
  std::vector<float> output;
  {
    ds::StagedVideoRequest<Filter> request(n,s);
    while(!request.advance({&s.source,1},s)) {
      CHECK(request.pending().size()==1);
      CHECK(live==0); // Previous acquired owner was released before the next stage.
      const auto index=request.pending()[0].frame_number;trace.sources.push_back(index);
      if(index==fail_source) throw std::runtime_error("injected delivery failure");
      auto pixels=std::make_shared<Pixels>(&live);
      for(int y=0;y<24;++y) for(int x=0;x<32;++x)
        pixels->pixels[y*32+x]=index==0 ? std::numeric_limits<float>::quiet_NaN() : float((std::int64_t(index)*7+y*3+x*11)%127)/255;
      ds::FrameRef owner(pixels);request.accept({0,index,owner.view(),owner});
      trace.peak=std::max(trace.peak,live);CHECK(live==1);
    }
    Factory factory;factory.fail=fail_output;
    auto frame=request.finish({32,24,s.source.num_frames,s.source.format,{}},{&s.source,1},s,factory);
    const auto view=frame.view();const auto* data=static_cast<const float*>(view.plane(0).data);
    output.assign(data,data+32*24);
  }
  CHECK(live==0);return output;
}
int main() {try {
  using Z=std::complex<float>;
  for(std::size_t count:{1u,2u,3u,7u,17u,33u,65u}) {
    std::vector<Z> x(count+2,{1,0}),l(count+2);
    std::vector<float> c(count+2,2),q=c;
    l.front()=l.back()={99,99};
    kalman_scalar(x.data()+1,l.data()+1,c.data()+1,q.data()+1,nullptr,2,4,count);
    CHECK(l.front()==Z(99,99) && l.back()==Z(99,99));
    for(std::size_t i=1;i<=count;++i) CHECK(l[i]==Z(4.f/6.f,0));
  }
  {Z x{2,0},l{};float c=2,q=2;
    kalman_scalar(&x,&l,&c,&q,nullptr,1,4,1);CHECK(l.real()!=2); // equality smooths
    x={0,3};l={0,0};c=q=2;kalman_scalar(&x,&l,&c,&q,nullptr,1,4,1);
    CHECK(l==x && c==1 && q==c); // imaginary-only motion resets both
    x={0,0};l={0,0};c=q=0;kalman_scalar(&x,&l,&c,&q,nullptr,0,4,1);CHECK(l==x && c==0);
    x={1e30f,0};l={-1e30f,0};c=q=2;kalman_scalar(&x,&l,&c,&q,nullptr,2,4,1);CHECK(l==x);
    x={0,0};l={0,0};c=q=2e38f;rejects([&]{kalman_scalar(&x,&l,&c,&q,nullptr,1,4,1);});
    x={NAN,0};rejects([&]{kalman_scalar(&x,&l,&c,&q,nullptr,0,4,1);});
  }
  {
    FFT3DConfig c;c.bt=0;c.bw=c.bh=8;c.opt=1;c.sigma=12;c.pfactor=.7f;c.px=c.py=2;
    Plan patterned(32,24,{8,false,false},c);auto initial=patterned.initial_kalman();
    CHECK(initial.last[0]==Z() && initial.covariance[0]==12.f*12.f*64.f);
    c.pfactor=0;c.enhancement.sharpen=.5f;Plan enhanced(32,24,{32,true,false},c);
    c.enhancement.sharpen=0;Plan plain(32,24,{32,true,false},c);
    auto a=enhanced.initial_kalman(),b=plain.initial_kalman();std::vector<float> samples(32*24,.25f),out(samples.size());
    span2d::Plane<const float> source(samples.data(),32,24,32*4);
    enhanced.advance_kalman(source,a);plain.advance_kalman(source,b);
    enhanced.render_kalman(source,span2d::Plane<float>(out.data(),32,24,32*4),a);
    CHECK(a.last==b.last && a.covariance==b.covariance && a.process==b.process);
    enhanced.advance_kalman(source,a);plain.advance_kalman(source,b);CHECK(a.last==b.last && a.covariance==b.covariance);
  }
  for(bool sampled:{false,true}) {
    FFT3DConfig config;config.sigma=12;config.pfactor=sampled ? 1 : 0;config.px=2;config.py=2;
    auto canonical=state(config,0);std::array<std::vector<float>,13> expected;
    for(int n=0;n<13;++n) {Trace trace;expected[n]=run(canonical,n,trace);CHECK(trace.peak==1);}
    auto cached=state(config);
    cached.checkpoints=std::make_shared<runtime::Checkpoints>();
    for(int n:{8,2,12,4,8,6,1,11,3,12}) {Trace trace;auto result=run(cached,n,trace);CHECK(result==expected[n]);}
    auto concurrent=state(config);
    std::vector<std::future<std::vector<float>>> jobs;
    for(int n:{12,3,12,8}) jobs.push_back(std::async(std::launch::async,[&,n]{Trace trace;return run(concurrent,n,trace);}));
    int i=0;for(int n:{12,3,12,8}) CHECK(jobs[i++].get()==expected[n]);
    auto failing=state(config);Trace trace;
    rejects([&]{run(failing,8,trace,4);});CHECK(!failing.checkpoints->acquire(8));
    CHECK(run(failing,8,trace)==expected[8]);
    auto allocation=state(config);rejects([&]{run(allocation,8,trace,-1,true);});CHECK(!allocation.checkpoints->acquire(8));
    CHECK(run(allocation,8,trace)==expected[8]);
    auto exact=cached.checkpoints->acquire(12);CHECK(exact && exact->frame==12);
    Trace warm;CHECK(run(cached,12,warm)==expected[12]);CHECK(warm.sources==std::vector<int>{12});
  }
  { // INT_MAX clip boundary with a synthetic immediately preceding checkpoint.
    auto s=state({});auto cp=std::make_unique<runtime::Checkpoint>();cp->frame=INT32_MAX-2;cp->planes[0]=s.plans[0]->initial_kalman();
    s.checkpoints->publish(std::move(cp));Trace trace;run(s,INT32_MAX-1,trace);CHECK(trace.sources==std::vector<int>{INT32_MAX-1});
  }
  {runtime::Checkpoints cache;auto first=std::make_unique<runtime::Checkpoint>();first->frame=1;cache.publish(std::move(first));
    auto lease=cache.acquire(1);
    for(int n=2;n<=6;++n){auto cp=std::make_unique<runtime::Checkpoint>();cp->frame=n;cache.publish(std::move(cp));}
    CHECK(lease->frame==1);CHECK(!cache.acquire(1));CHECK(cache.bytes()<=64*1024*1024);
  }
  { // Cache admission uses actual capacity, including ownership/cache metadata.
    constexpr std::size_t budget=64*1024*1024;
    auto large=[](int frame) {
      auto cp=std::make_unique<runtime::Checkpoint>();cp->frame=frame;
      cp->planes[0].last.reserve(budget/sizeof(Z));
      return cp;
    };
    runtime::Checkpoints cache;
    cache.publish(large(1));auto lease=cache.acquire(1);
    CHECK(lease && cache.bytes()==sizeof(cache)+lease->bytes());
    CHECK(cache.bytes()>budget);
    cache.publish(large(2));auto next=cache.acquire(2);
    CHECK(next && next->frame==2 && !cache.acquire(1));
    CHECK(cache.bytes()==sizeof(cache)+next->bytes()); // Exactly one oversized state.
    CHECK(lease->frame==1); // Eviction cannot invalidate an in-flight request.
    runtime::Checkpoints hard(budget);hard.publish(large(1));CHECK(!hard.acquire(1));
    runtime::Checkpoints zero(0);auto cp=std::make_unique<runtime::Checkpoint>();cp->frame=1;
    zero.publish(std::move(cp));CHECK(!zero.acquire(1));
  }
  { // Account for mixed element sizes and reserved (not only live) storage.
    KalmanState s;s.last.reserve(11);s.covariance.reserve(17);s.process.reserve(19);
    CHECK(s.bytes()==sizeof(s)+s.last.capacity()*sizeof(Z)+
        (s.covariance.capacity()+s.process.capacity())*sizeof(float));
    FFT3DConfig c;c.bt=0;c.bw=c.bh=32;c.ow=c.oh=16;c.opt=1;
    auto cp=std::make_unique<runtime::Checkpoint>();cp->frame=1;std::size_t payload=0;
    for(int p=0;p<3;++p) {
      Plan plan(p ? 640:1280,p ? 360:720,{8,false,false},c);
      cp->planes[p]=plan.initial_kalman();
      payload+=cp->planes[p].last.size()*16;
    }
    CHECK(payload==49560576); // 47.265 MiB for 720p YUV420, formerly 70.897 MiB.
    runtime::Checkpoints hard(64*1024*1024);hard.publish(std::move(cp));
    CHECK(hard.acquire(1) && hard.bytes()<=64*1024*1024);
  }
  std::cout<<"Kalman scalar, canonical replay/cache/failures, exact checkpoint, source-owner bound and INT_MAX passed\n";
  return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
