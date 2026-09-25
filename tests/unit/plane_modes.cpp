#include "plugin/filter.hpp"
#include "../test.hpp"
using namespace neo_fft;

// Defined sentinel storage lets us test no writes without reading uninitialized
// host memory. Padding also detects writes outside each selected plane.
template<class T> struct Pixels : ds::FrameStorage {
  ds::VideoFormat format;
  std::array<std::vector<T>,4> data;
  explicit Pixels(ds::VideoFormat f, int frame=-1) : format(f) {
    for(int p=0;p<f.plane_count;++p) {
      const int w=width(p),h=height(p);
      data[p].assign((w+4)*h,T(201));
      if(frame>=0) for(int y=0;y<h;++y) for(int x=0;x<w;++x) {
        const int value=16+(frame*17+p*23+y*7+x*11)%180;
        data[p][y*(w+4)+x]=std::is_same_v<T,float> ? T(value/255.f) : T(value);
      }
    }
  }
  bool chroma(int p)const{return format.color_family==ds::ColorFamily::Yuv && (p==1 || p==2);}
  int width(int p)const{return 32>>(chroma(p) ? format.subsampling_w : 0);}
  int height(int p)const{return 24>>(chroma(p) ? format.subsampling_h : 0);}
  ds::VideoFrameView read()const override {
    ds::VideoFrameView v;v.format=format;v.plane_count=format.plane_count;
    for(int p=0;p<v.plane_count;++p)v.planes[p]={data[p].data(),std::ptrdiff_t((width(p)+4)*sizeof(T)),width(p),height(p)};
    return v;
  }
  ds::MutableVideoFrameView write()override {
    ds::MutableVideoFrameView v;v.format=format;v.plane_count=format.plane_count;
    for(int p=0;p<v.plane_count;++p)v.planes[p]={data[p].data(),std::ptrdiff_t((width(p)+4)*sizeof(T)),width(p),height(p)};
    return v;
  }
};
template<class T> struct Factory : ds::FrameFactory {
  ds::WritableFrame allocate(ds::VideoFormat f,int,int,const ds::FrameRef&)override {
    return ds::WritableFrame(std::make_unique<Pixels<T>>(f));
  }
};
template<Algorithm A,class T>
ds::FrameRef render(typename plugin::Filter<A>::State& state,int n) {
  using F=plugin::Filter<A>;
  ds::StagedVideoRequest<F> request(n,state);
  while(!request.advance({&state.source,1},state)) {
    for(const auto& pending:request.pending()) {
      ds::FrameRef owner(std::make_shared<Pixels<T>>(state.source.format,pending.frame_number));
      request.accept({0,pending.frame_number,owner.view(),owner});
    }
  }
  Factory<T> factory;
  return request.finish({32,24,8,state.source.format,{}},{&state.source,1},state,factory);
}
template<Algorithm A,class T>
void check(ds::VideoFormat format,ds::ParamValues params) {
  using F=plugin::Filter<A>;
  ds::VideoInputInfo input{32,24,8,format};
  constexpr const char* names[]={"y","u","v","a"};
  for(int rotation=0;rotation<5;++rotation) {
    auto mixed=params;
    std::vector<int> processed;
    std::array<int,4> modes{};
    for(int p=0;p<4;++p) {
      modes[p]=rotation==4 ? 1 : 1+(p+rotation)%3;
      mixed.entries.push_back({names[p],modes[p]});
      if(p<format.plane_count && modes[p]==3)processed.push_back(p);
    }
    ds::VideoInitContext ctx;ctx.inputs={&input,1};ctx.params=&mixed;
    auto actual=plugin::unwrap(F::init(ctx)).state;
    auto reference_params=mixed;
    reference_params.entries.push_back({"planes",processed});
    ctx.params=&reference_params;
    auto reference=plugin::unwrap(F::init(ctx)).state;
    for(int n:{0,4,2,4}) {
      const auto output=render<A,T>(actual,n),expected=render<A,T>(reference,n);
      const auto view=output.view(),oracle=expected.view();
      for(int p=0;p<format.plane_count;++p) {
        const auto& plane=view.plane(p);const auto& ref=oracle.plane(p);
        const auto* data=static_cast<const T*>(plane.data);
        const auto* wanted=static_cast<const T*>(ref.data);
        for(int y=0;y<plane.height;++y)for(int x=0;x<plane.width+4;++x) {
          const auto offset=y*(plane.width+4)+x;
          if(modes[p]==1 || x>=plane.width)CHECK(data[offset]==T(201));
          else CHECK(data[offset]==wanted[offset]);
        }
      }
    }
  }
}
template<class T> void cases(ds::SampleFormat sample) {
  for(int opt:{0,1})for(auto family:{ds::ColorFamily::Gray,ds::ColorFamily::Yuv,ds::ColorFamily::Rgb}) {
    const bool yuv=family==ds::ColorFamily::Yuv;
    const ds::VideoFormat f{family,sample,family==ds::ColorFamily::Gray ? 1 : 4,yuv ? 1 : 0,yuv ? 1 : 0};
    for(int bt:{0,3}) {
      ds::ParamValues p;p.entries={{"bw",8},{"bh",8},{"bt",bt},{"opt",opt}};
      check<Algorithm::FFT3D,T>(f,p);
      p.entries.insert(p.entries.end(),{{"interlaced",true},{"l",2},{"r",2}});
      check<Algorithm::FFT3D,T>(f,p);
    }
    for(int temporal:{0,1}) {
      ds::ParamValues p;p.entries={{"sbsize",4},{"sosize",2},{"tbsize",temporal ? 4 : 3},
        {"tmode",temporal},{"tosize",temporal ? 2 : 0},{"opt",opt}};
      check<Algorithm::DFTTest,T>(f,p);
    }
  }
}
int main(){try {
  cases<std::uint8_t>(ds::SampleFormat::UInt8);
  cases<std::uint16_t>(ds::SampleFormat::UInt16);
  cases<float>(ds::SampleFormat::Float32);
  std::cout<<"Plane skip sentinels, mixed modes, precedence, ROI and Kalman replay passed\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
