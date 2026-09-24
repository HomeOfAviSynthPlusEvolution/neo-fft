#include "plugin/filter.hpp"
#include "../test.hpp"
using namespace neo_fft;
using Filter=plugin::Filter<Algorithm::DFTTest>;
struct Pixels:ds::FrameStorage {
  std::vector<float> data=std::vector<float>(32*24);int* live;
  explicit Pixels(int* count=nullptr):live(count){if(live)++*live;}
  ~Pixels()override{if(live)--*live;}
  ds::VideoFrameView read()const override{return {{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0},1,{{{data.data(),128,32,24}}}};}
  ds::MutableVideoFrameView write()override{return {{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0},1,{{{data.data(),128,32,24}}}};}
};
struct Factory:ds::FrameFactory {
  bool fail=false;
  ds::WritableFrame allocate(ds::VideoFormat,int,int,const ds::FrameRef&)override {
    if(fail)throw std::bad_alloc();return ds::WritableFrame(std::make_unique<Pixels>());
  }
};
Filter::State make_state(bool sampled=false) {
  DFTConfig c;c.temporal_mode=1;c.tbsize=4;c.temporal_overlap=2;c.block=4;c.overlap=2;c.opt=1;
  if(sampled)c.locations={{1,0,0,0}};
  Filter::State state;state.source={32,24,INT32_MAX,{ds::ColorFamily::Gray,ds::SampleFormat::Float32,1,0,0}};
  state.temporal_mode=1;state.temporal_size=4;state.temporal_overlap=2;
  state.plans[0]=std::make_shared<Plan>(32,24,SampleFormat{32,true,false},c);state.dft_noise=state.plans[0]->dft_noise();state.sample_bits=32;
  return state;
}
std::vector<float> run(Filter::State& state,int n,int& live,int fail=-1,bool fail_output=false,bool cancel=false) {
  ds::StagedVideoRequest<Filter> request(n,state);
  CHECK(!request.advance({&state.source,1},state));
  CHECK(request.pending().size()<=std::size_t(7+(state.dft_noise ? 4 : 0)));
  std::vector<int> acquired;
  for(const auto& pending:request.pending()) {
    const int index=pending.frame_number;CHECK(std::find(acquired.begin(),acquired.end(),index)==acquired.end());acquired.push_back(index);
    auto pixels=std::make_shared<Pixels>(&live);
    for(int y=0;y<24;++y)for(int x=0;x<32;++x)pixels->data[y*32+x]=float((std::int64_t(index)*17+y*7+x*11)%127)/255;
    if(index==fail)pixels->data[0]=NAN;
    ds::FrameRef owner(pixels);request.accept({0,index,owner.view(),owner});
  }
  std::sort(acquired.begin(),acquired.end());
  std::vector<int> expected;
  for(std::int64_t index=std::int64_t(n)-2;index<=std::min(std::int64_t(n)+3,std::int64_t(INT32_MAX)-1);++index)expected.push_back(int(index));
  if(state.dft_noise)expected.insert(expected.end(),{1,2,3,4});
  std::sort(expected.begin(),expected.end());
  CHECK(acquired==expected);
  if(cancel)return {};
  CHECK(request.advance({&state.source,1},state));
  Factory factory;factory.fail=fail_output;
  auto out=request.finish({32,24,state.source.num_frames,state.source.format,{}},{&state.source,1},state,factory);
  auto view=out.view();const auto* p=static_cast<const float*>(view.plane(0).data);return {p,p+32*24};
}
void unique_fetches() {
  struct Provider:ds::VideoFrameProvider {
    Pixels source;
    std::vector<int> fetched;
    ds::Result<ds::RequestedVideoFrame> get(int,int n)override {
      fetched.push_back(n);
      return ds::Result<ds::RequestedVideoFrame>::success({0,n,source.read(),{}});
    }
  } provider;
  auto state=make_state();
  state.source.num_frames=100;
  Pixels output;
  for(const auto& [n,expected]:std::vector<std::pair<int,std::vector<int>>>{{0,{0,1,2,3}},{8,{6,7,8,9,10,11}}}) {
    provider.fetched.clear();
    ds::VideoProcessContext ctx{n,provider,output.write(),&state};
    CHECK(Filter::process(ctx).has_value());
    CHECK(provider.fetched==expected);
  }
  DFTConfig centered;centered.tbsize=5;centered.block=4;centered.overlap=2;centered.opt=1;
  state.temporal_mode=0;state.temporal_size=5;
  state.plans[0]=std::make_shared<Plan>(32,24,SampleFormat{32,true,false},centered);
  for(const auto& [n,expected]:std::vector<std::pair<int,std::vector<int>>>{{0,{0,1,2}},{8,{6,7,8,9,10}}}) {
    provider.fetched.clear();
    ds::VideoProcessContext ctx{n,provider,output.write(),&state};
    CHECK(Filter::process(ctx).has_value());
    CHECK(provider.fetched==expected);
  }
}
void unique_requests() {
  auto state=make_state(true);
  state.source.num_frames=100;
  const auto check=[&](int n,std::vector<ds::VideoFrameRequest> requests,const std::vector<int>& expected) {
    ds::VideoRequestContext ctx{n,requests,{},&state};
    CHECK(Filter::request(ctx).has_value());
    std::vector<int> actual;
    for(const auto& request:requests) {CHECK(request.input_index==0);actual.push_back(request.frame_number);}
    CHECK(actual==expected);
  };
  check(8,{}, {1,2,3,4,6,7,8,9,10,11});
  check(8,{{0,8}}, {8,1,2,3,4,6,7,9,10,11});
  DFTConfig centered;centered.tbsize=5;centered.block=4;centered.overlap=2;centered.opt=1;
  state.temporal_mode=0;state.temporal_size=5;state.dft_noise.reset();
  state.plans[0]=std::make_shared<Plan>(32,24,SampleFormat{32,true,false},centered);
  check(0,{}, {0,1,2});
}
template<Algorithm A> void plane_admission() {
  using F=plugin::Filter<A>;
  ds::VideoInputInfo input{128,96,8,{ds::ColorFamily::Yuv,ds::SampleFormat::UInt8,4,1,1}};
  ds::ParamValues params;
  ds::VideoInitContext ctx;ctx.inputs={&input,1};ctx.params=&params;
  auto ordinary=plugin::unwrap(F::init(ctx)).state;
  CHECK(ordinary.plans[0] && ordinary.plans[1] && ordinary.plans[2] && !ordinary.plans[3]);
  params.entries={{"y",2},{"u",2},{"v",2},{"a",2}};
  // All-copy must not enter plan geometry/model setup, even with tiny planes.
  input.width=input.height=4;
  auto copied=plugin::unwrap(F::init(ctx)).state;
  for(const auto& plan:copied.plans)CHECK(!plan);
  CHECK(!copied.dft_noise && !copied.kalman && !copied.sampled);
  std::vector<ds::VideoFrameRequest> requests;
  ds::VideoRequestContext request{4,requests,{},&copied};
  CHECK(F::request(request).has_value());
  CHECK(requests.size()==1 && requests[0].frame_number==4);
  input.width=128;input.height=96;
  params.entries.push_back({"planes",std::vector<int>{3}});
  auto alpha=plugin::unwrap(F::init(ctx)).state;
  CHECK(alpha.plans[3] && !alpha.plans[0] && !alpha.plans[1] && !alpha.plans[2]);
}
int main(){try {
  plane_admission<Algorithm::FFT3D>();plane_admission<Algorithm::DFTTest>();
  unique_fetches();
  unique_requests();
  for(bool sample:{false,true}) {
    int live=0;auto state=make_state(sample),clean=make_state(sample);
    const auto expected=run(clean,8,live);CHECK(live==0);
    rejects([&]{run(state,8,live,sample ? 1 : 8);});CHECK(live==0);
    if(sample)CHECK(!state.dft_noise->power());
    CHECK(run(state,8,live)==expected && live==0);
    rejects([&]{run(state,8,live,-1,true);});CHECK(live==0);
    CHECK(run(state,8,live)==expected && live==0);
    run(state,8,live,-1,false,true);CHECK(live==0);
    const auto bytes=state.plans[0]->workspace_pool().budget().total_bytes;
    run(state,INT32_MAX-1,live);CHECK(live==0);
    CHECK(state.plans[0]->workspace_pool().budget().total_bytes==bytes);
  }
  std::cout<<"Temporal closure/owners, cancellation, model/output failure retry and far-index bounds passed\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
