#include "plugin/filter.hpp"
#include "../test.hpp"
#include <future>
#include <thread>
using namespace neo_fft;
using Cache=runtime::SpectraCache;

void ownership() {
  const std::array<std::size_t,4> bins{64000,32000,32000};
  Cache cache(bins,3,100,1,{3,3,3});
  const auto fill=[](std::complex<float>* dst){dst[0]={.25f,.5f};};
  auto a=cache.get(0,0,0,fill),b=cache.get(0,0,1,fill);
  CHECK(a && b && a.get()!=b.get() && cache.bytes()>1024000 && cache.bytes()<=1048576);
  CHECK(!cache.get(1,0,0,fill)); // Both rows are pinned.
  a.reset();auto other_plane=cache.get(0,1,0,fill);CHECK(other_plane);
  CHECK(!cache.get(1,0,0,fill));
  other_plane.reset();auto c=cache.get(1,0,0,fill);CHECK(c);
  CHECK(b.get()[0]==std::complex<float>(.25f,.5f));CHECK(cache.peak_bytes()<=1048576);
  CHECK(cache.get(0,0,1,fill).get()==b.get());
  rejects([&]{cache.get(0,0,3,fill);});
  const auto stats=cache.stats();
  CHECK(stats.requests==7 && stats.builds==4 && stats.bypasses==2 && stats.hits==1);
  CHECK(stats.frames==2 && stats.rows==2 && stats.evictions==2);
  Cache small(bins,3,1,1);
  small.inform_threads(4);CHECK(small.frame_limit()==1);
  auto one=small.get(0,0,0,fill);CHECK(one && !small.get(1,0,0,fill));
  one.reset();CHECK(small.get(1,0,0,fill));
  Cache row_frames({8,8,0},3,1,1,{3,2,0});
  auto row_pin=row_frames.get(0,0,0,fill);
  CHECK(row_frames.get(0,0,1,fill) && row_frames.get(0,1,0,fill));
  CHECK(row_frames.stats().rows==3 && row_frames.stats().frames==1);
  CHECK(!row_frames.get(1,0,0,fill)); // One pinned row protects its frame from frame-limit eviction.
  row_pin.reset();CHECK(row_frames.get(1,0,0,fill));
  CHECK(row_frames.stats().rows==1 && row_frames.stats().frames==1 && row_frames.stats().evictions==3);
  Cache off(bins,3,0,512),off_mb(bins,3,9,0),spatial(bins,1),large({200000,0,0},3,9,1);
  CHECK(!off.get(0,0,0,fill) && !off_mb.get(0,0,0,fill) && !spatial.get(0,0,0,fill) && !large.get(0,0,0,fill));
  CHECK(large.bytes()==0);
  rejects([&]{Cache invalid(bins,3,-2);});rejects([&]{Cache invalid(bins,3,1,-2);});
  using F=plugin::Filter<Algorithm::FFT3D>;
  F::State state;state.spectra=std::make_shared<Cache>(bins,3);
  CHECK(state.spectra->frame_limit()==3);
  ds::VideoCacheHintsContext hint{514,4,0,&state};F::cache_hints(hint);CHECK(state.spectra->frame_limit()==6);
  hint.frame_range=1;F::cache_hints(hint);CHECK(state.spectra->frame_limit()==3);
  hint.frame_range=0;F::cache_hints(hint);CHECK(state.spectra->frame_limit()==3);
  Cache resized({8,0,0},3);
  resized.inform_threads(4);
  std::array<Cache::Lease,6> pins;
  for(int n=0;n<6;++n)pins[n]=resized.get(n,0,0,fill);
  const auto before=resized.bytes();resized.inform_threads(1);
  CHECK(resized.bytes()==before);pins={};
  CHECK(resized.get(5,0,0,fill));CHECK(resized.bytes()==before/2);
}

void alpha_rows() {
  Cache cache({8,0,0,8},3,3,1,{2,0,0,2});
  auto registration=cache.register_request(0,2,{false,false,false,true});
  const auto luma=cache.get(1,0,0,[](auto* data){data[0]={1,2};});
  const auto alpha=cache.get(1,3,0,[](auto* data){data[0]={3,4};});
  CHECK(luma && alpha && luma.get()!=alpha.get());
  CHECK(cache.get(1,3,0,[](auto*){throw std::runtime_error("Alpha cache miss");}).get()==alpha.get());
  registration->complete_row(3,0);
  CHECK(cache.get(1,3,1,[](auto* data){data[0]={5,6};}));
  CHECK(cache.stats().builds==3 && cache.stats().hits==1);
  rejects([&]{cache.get(1,4,0,[](auto*){});});
  rejects([&]{cache.get(1,3,2,[](auto*){});});
}

void registered_eviction() {
  const auto fill=[](auto* dst){dst[0]={7,9};};
  const auto must_hit=[](auto*){throw std::runtime_error("evicted pending row");};
  Cache cache({64000,0,0},3,100,1,{4,0,0});
  auto request=cache.register_request(0,0);
  CHECK(cache.get(0,0,0,fill) && cache.get(1,0,0,fill));
  CHECK(cache.get(2,0,0,fill)); // Keep old frame 0 for its registered consumer.
  CHECK(cache.get(0,0,0,must_hit));
  CHECK(cache.stats().evicted_needed_rows==0 && cache.stats().active_requests==1);
  request->complete_row(0,0);
  CHECK(cache.get(3,0,0,fill)); // Frame 0 is now expendable even though just touched.
  int rebuilds=0;
  CHECK(cache.get(0,0,0,[&](auto* out){++rebuilds;fill(out);}));CHECK(rebuilds==1);
  request.reset();CHECK(cache.stats().active_requests==0);

  Cache rows({64000,0,0},3,100,1,{3,0,0});
  auto progress=rows.register_request(0,0);
  CHECK(rows.get(0,0,1,fill) && rows.get(0,0,0,fill));
  progress->complete_row(0,0);
  CHECK(rows.get(1,0,0,fill) && rows.get(0,0,1,must_hit));
  CHECK(rows.stats().evicted_needed_rows==0);

  Cache pressure({64000,0,0},3,100,1);
  auto all=pressure.register_request(0,2);
  CHECK(pressure.get(0,0,0,fill));CHECK(pressure.get(1,0,0,fill));CHECK(pressure.get(2,0,0,fill));
  CHECK(pressure.stats().evicted_needed_rows==1 && pressure.stats().evicted_consumers==1);
  CHECK(pressure.peak_bytes()<=1048576); // Future demand is not an unbounded pin.
  auto duplicate=pressure.register_request(0,2);
  CHECK(pressure.get(3,0,0,fill));
  CHECK(pressure.stats().evicted_needed_rows==2 && pressure.stats().evicted_consumers==3);
  duplicate.reset();all.reset();CHECK(pressure.stats().active_requests==0);
  try {auto cancelled=pressure.register_request(9,9);throw std::runtime_error("cancel");}
  catch(const std::runtime_error&) {}
  CHECK(pressure.stats().active_requests==0 && pressure.stats().peak_requests==2);
  rejects([&]{pressure.register_request(-1,0);});rejects([&]{pressure.register_request(0,3);});

  Cache frames({8,8,0},3,2,1,{2,2,0});
  auto pending=frames.register_request(0,0);
  CHECK(frames.get(0,0,0,fill) && frames.get(0,1,1,fill) && frames.get(1,0,0,fill));
  CHECK(frames.get(2,0,0,fill));
  CHECK(frames.get(0,0,0,must_hit) && frames.get(0,1,1,must_hit));
  pending->complete_row(0,1);pending->complete_row(1,1);
  CHECK(frames.get(3,0,0,fill));
  CHECK(frames.stats().frames==2 && frames.stats().rows==2 && frames.stats().evicted_needed_rows==0);
}

void staged_registration() {
  using F=plugin::Filter<Algorithm::FFT3D>;
  for(int bt:{2,3,4,5}) {
    F::State state;state.source={32,24,11,{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0}};
    FFT3DConfig config;config.bt=bt;config.bw=config.bh=8;config.ow=config.oh=4;config.opt=1;
    state.temporal_size=bt;state.plans[0]=std::make_shared<Plan>(32,24,SampleFormat{32,true,false},config);
    const auto& plan=*state.plans[0];
    state.spectra=std::make_shared<Cache>(std::array<std::size_t,4>{std::size_t(plan.geometry.x.count)*plan.fft.bins(),0,0},
        bt,10,1,std::array<int,4>{plan.geometry.y.count,0,0});
    for(int n:{0,5,10}) {
      {
        ds::StagedVideoRequest<F> request(n,state);
        CHECK(!request.advance({&state.source,1},state));
        const bool edge=n<bt/2 || 10-n<(bt-1)/2;
        CHECK(request.pending().size()==std::size_t(edge ? 1 : bt));
        CHECK(state.spectra->stats().active_requests==1); // Registered before any source delivery.
        {
          ds::StagedVideoRequest<F> duplicate(n,state);
          CHECK(!duplicate.advance({&state.source,1},state));
          CHECK(state.spectra->stats().active_requests==2);
        }
        CHECK(state.spectra->stats().active_requests==1);
      } // Host cancellation before dependency delivery removes its registration.
      CHECK(state.spectra->stats().active_requests==0);
    }
  }
}

void publication(bool fail) {
  Cache cache({8,0,0},3,1,1);
  std::promise<void> entered,release;auto gate=release.get_future().share();
  auto produce=std::async(std::launch::async,[&]{
    return cache.get(0,0,0,[&](auto* out){entered.set_value();CHECK(gate.wait_for(std::chrono::seconds(5))==std::future_status::ready);if(fail)throw std::runtime_error("producer failure");out[0]={7,9};});
  });
  entered.get_future().wait();
  auto consume=[&]{return cache.get(0,0,0,[](auto*){throw std::runtime_error("duplicate producer");});};
  auto b=std::async(std::launch::async,consume),c=std::async(std::launch::async,consume);
  const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
  while(cache.waits()<2 && std::chrono::steady_clock::now()<deadline)std::this_thread::yield();
  const bool both_waiting=cache.waits()==2;
  bool fallback=false;
  try {fallback=!cache.get(1,0,0,[](auto*){});} // Never waits for budget held by a builder.
  catch(...) {release.set_value();throw;}
  release.set_value();CHECK(both_waiting && fallback);
  if(fail) {
    auto producer_error=[](auto& future) {
      bool caught=false;
      try {future.get();}catch(const std::runtime_error& e) {caught=std::string(e.what())=="producer failure";}
      CHECK(caught);
    };
    producer_error(produce);producer_error(b);producer_error(c);
    CHECK(!cache.get(0,0,0,[](auto*){}));
    CHECK(cache.get(1,0,0,[](auto* p){p[0]={1,0};})); // Evict failure, retry cleanly.
    CHECK(cache.get(0,0,0,[](auto* p){p[0]={7,9};}));
  } else {
    auto x=produce.get(),y=b.get(),z=c.get();
    CHECK(x.get()==y.get() && y.get()==z.get() && x.get()[0]==std::complex<float>(7,9));
    CHECK(cache.builds()==1 && cache.hits()==2);
  }
}

template<class T> void reconstruction(SampleFormat format) {
  constexpr int w=37,h=29,count=11;
  constexpr std::ptrdiff_t stride=w*static_cast<std::ptrdiff_t>(sizeof(T));
  std::array<std::vector<T>,count> input;
  for(int n=0;n<count;++n) {
    input[n].resize(w*h);
    for(int k=0;k<w*h;++k) {
      if constexpr(std::is_same_v<T,float>)input[n][k]=float((k*31+n*17)%251-125)/255;
      else input[n][k]=T(((k*31+n*17)%251)*(1<<(format.bits-8)));
    }
  }
  for(int opt:{0,1})for(int bt:{1,2,3,4,5}) {
    FFT3DConfig config;config.bw=config.bh=8;config.ow=config.oh=4;config.bt=bt;config.opt=opt;
    config.sigma2=4;config.enhancement.sharpen=.3f;
    Plan plan(w,h,format,config);
    const auto bins=std::size_t(plan.geometry.x.count)*plan.fft.bins();
    auto render=[&](int n,Cache* cache) {
      const bool edge=n<bt/2 || count-1-n<(bt-1)/2;
      const int slots=edge ? 1 : bt;
      std::vector<span2d::Plane<const T>> sources;
      for(int j=0;j<slots;++j)sources.emplace_back(input[n-slots/2+j].data(),w,h,stride);
      std::vector<T> result(w*h);
      auto registration=cache ? cache->register_request(n-slots/2,n-slots/2+slots-1,{true,false,false}) : nullptr;
      plan.process_at<T>({sources.data(),sources.size()},{result.data(),w,h,stride},n,0,{},cache,registration.get());
      return result;
    };
    std::array<std::vector<T>,count> expected;
    for(int n=0;n<count;++n)expected[n]=render(n,nullptr);
    for(int frames:{0,1,3,8}) {
      Cache cache({bins,0,0},bt,frames,1,{plan.geometry.y.count,0,0});
      for(int n:{0,1,2,3,4,5,6,7,8,9,10,9,4,1,7,3,8,6,5,4,0}) {
        auto actual=render(n,&cache);
        CHECK(std::memcmp(actual.data(),expected[n].data(),actual.size()*sizeof(T))==0);
      }
      for(int n:{3,5,8}) {
        auto a=std::async(std::launch::async,[&]{return render(n,&cache);});
        auto b=std::async(std::launch::async,[&]{return render(n,&cache);});
        auto c=std::async(std::launch::async,[&]{return render(n-1,&cache);});
        CHECK(a.get()==expected[n] && b.get()==expected[n] && c.get()==expected[n-1]);
      }
      CHECK(cache.peak_bytes()<=1048576 && cache.stats().active_requests==0);
      if(bt>1 && frames>=3)CHECK(cache.hits()>0 && cache.builds()>0);
    }
    if constexpr(std::is_same_v<T,float>)if(bt==3) {
      Cache cache({bins,0,0},bt,8,1,{plan.geometry.y.count,0,0});
      const auto saved=input[4][0];input[4][0]=NAN;
      rejects([&]{render(4,&cache);});
      input[4][0]=saved;
      CHECK(render(4,&cache)==expected[4]);
      CHECK(cache.peak_bytes()<=1048576 && cache.stats().active_requests==0);
    }
  }
}
int main() {try {
  ownership();alpha_rows();registered_eviction();staged_registration();publication(false);publication(true);
  reconstruction<std::uint8_t>({8,false,false});
  reconstruction<std::uint16_t>({16,false,true});
  reconstruction<float>({32,true,true});
  std::cout<<"spectra cache: budgets, hints, publication, failures and ordered/concurrent reconstruction passed\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
