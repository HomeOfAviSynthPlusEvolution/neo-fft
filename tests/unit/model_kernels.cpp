#include "kernels/model.hpp"
#include "kernels/table.hpp"
#include "kernels/spatial.hpp"
#include "algorithms/windows.hpp"
#include "base/checked.hpp"
#include "../test.hpp"
#include <cstring>
#ifdef NEO_FFT_TEST_HIGHWAY
#include <hwy/targets.h>
#endif
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
using namespace neo_fft;
template<class T> struct Guard {
  T* data;
  void* allocation=nullptr;
  std::vector<T> fallback;
  explicit Guard(std::size_t n) {
#if defined(_WIN32)
    SYSTEM_INFO info{};GetSystemInfo(&info);CHECK(n*sizeof(T)<=info.dwPageSize);
    allocation=VirtualAlloc(nullptr,2*info.dwPageSize,MEM_RESERVE|MEM_COMMIT,PAGE_READWRITE);CHECK(allocation);
    auto* boundary=static_cast<unsigned char*>(allocation)+info.dwPageSize;DWORD old;
    CHECK(VirtualProtect(boundary,info.dwPageSize,PAGE_NOACCESS,&old));
    data=reinterpret_cast<T*>(boundary-n*sizeof(T));
    for(std::size_t i=0;i<n;++i) new(data+i) T{};
#else
    fallback.resize(n);data=fallback.data();
#endif
  }
  ~Guard() {
#if defined(_WIN32)
    VirtualFree(allocation,0,MEM_RELEASE);
#endif
  }
  Guard(const Guard&)=delete;
};
template<class T> void check_decode(std::size_t n,const ModelKernels& native) {
  Guard<T> src(n);Guard<float> a(n),b(n);
  for(std::size_t i=0;i<n;++i) {
    if constexpr(std::is_same_v<T,float>) {
      const float values[]={-.5f,-0.f,0.f,1e-30f,123.5f};src.data[i]=values[i%5];
    } else {
      const unsigned values[]={65535,32768,32767,255,256,0};src.data[i]=T(values[i%6]);
    }
  }
  for(float base:{0.f,128.f}) for(float scale:{1.f,255.f,1.f/256}) {
    native.decode(src.data,sample_storage<T>,a.data,n,base,scale);
    model_scalar().decode(src.data,sample_storage<T>,b.data,n,base,scale);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
  }
}
void run_tests() {
  const auto native=select_model(0);const auto scalar=model_scalar();
  const auto tables=select_table(0),reference=table_scalar();
  // Independent input frequencies exercise all segments and exact knots, with
  // a guarded logical tail for row inputs and outputs.
  for(std::size_t n=1;n<=257;++n) {
    Guard<float> x(n),a(n),b(n),c(n),e(n);
    for(std::size_t i=0;i<n;++i)x.data[i]=float(i%33)/32;
    tables.radius(a.data,x.data,n,.3f,3);reference.radius(b.data,x.data,n,.3f,3);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    const float sigmas[]={2,3,4,5};
    tables.analytic(a.data,x.data,n,sigmas,.001f);reference.analytic(b.data,x.data,n,sigmas,.001f);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    const float knots[]={0,.125f,.5f,.875f,1},values[]={3,17,0,21,4};
    tables.curve(a.data,x.data,n,knots,values,5,.37f);reference.curve(b.data,x.data,n,knots,values,5,.37f);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    tables.product(a.data,x.data,n,.39f,.71f);reference.product(b.data,x.data,n,.39f,.71f);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    tables.weights(a.data,x.data,n,.3f,.01f);reference.weights(b.data,x.data,n,.3f,.01f);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    for(bool sharp:{false,true})for(bool halo:{false,true}) {
      tables.enhancement(sharp ? a.data : nullptr,halo ? c.data : nullptr,x.data,n,.1f,.18f,2);
      reference.enhancement(sharp ? b.data : nullptr,halo ? e.data : nullptr,x.data,n,.1f,.18f,2);
      if(sharp)CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
      if(halo)CHECK(std::memcmp(c.data,e.data,n*sizeof(float))==0);
    }
    tables.divide(a.data,n,.13f);reference.divide(b.data,n,.13f);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    Guard<double> axis(n);for(std::size_t i=0;i<n;++i)axis.data[i]=double(i%31)/31;
    CHECK(tables.window(a.data,axis.data,n,.39,.71,.125f)==reference.window(b.data,axis.data,n,.39,.71,.125f));
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    CHECK(tables.maximum(a.data,n)==reference.maximum(b.data,n));
    tables.check_scale(a.data,n,15);reference.check_scale(b.data,n,15);
  }
  for(int id=0;id<12;++id)for(int block:{16,31,32})for(int mode:{0,1}) {
    const auto a=dft_window_3d(3,block,block/2,mode,id,id,2.5f,2.5f,0);
    const auto b=dft_window_3d(3,block,block/2,mode,id,id,2.5f,2.5f,1);
    CHECK(a.h==b.h && a.wscale==b.wscale);
  }
  for(auto k:{tables,reference}) {
    float f[32]{},out[32];
    // Overflow in an unused interpolation branch must not reject a valid bin.
    const float large[]={1,1,3e38f,3e38f};
    std::fill_n(f,32,1.f);k.analytic(out,f,32,large,1);
    const float huge[]={1e30f,1e30f,1e30f,1e30f};
    rejects([&]{k.analytic(out,f,32,huge,1);});
    const float knots[]={0,.5f,1},values[]={0,1e30f,1};
    std::fill_n(f,32,.5f);k.curve(out,f,32,knots,values,3,1);
    f[31]=NAN;rejects([&]{k.curve(out,f,32,knots,values,3,1);});
    std::fill_n(f,32,1.f);rejects([&]{k.enhancement(nullptr,out,f,32,1,.18f,1e30f);});
    double axis[32];std::fill_n(axis,32,1.0);
    rejects([&]{k.window(out,axis,32,1e40,1,0);});
    CHECK(k.window(nullptr,nullptr,0,1,1,.25f)==.25f);
    CHECK(k.maximum(nullptr,0)==0);k.check_scale(nullptr,0,1);
    for(std::size_t i=0;i<32;++i) {
      std::fill_n(f,32,1.f);f[i]=1e30f;
      rejects([&]{k.check_scale(f,32,1e20f);});
    }
  }
  for(auto k:{native,scalar}) {
    k.decode(nullptr,SampleStorage::F32,nullptr,0,0,1);
    k.window(nullptr,nullptr,0,1);k.power(nullptr,nullptr,0,nullptr,0,false);
    CHECK(k.score(nullptr,nullptr,0)==0);k.scale(nullptr,nullptr,0,1,1);
  }
  for(std::size_t n=1;n<=257;++n) {
    check_decode<std::uint8_t>(n,native);check_decode<std::uint16_t>(n,native);check_decode<float>(n,native);
    Guard<float> a(n),b(n),weights(n);Guard<std::complex<float>> spectrum(n),grid(n);
    for(std::size_t i=0;i<n;++i) {
      weights.data[i]=float(i%17)*.125f;
      spectrum.data[i]={float(i%31)-12.3f,float(i%13)+.7f};
      grid.data[i]={float(i%7),float(i%5)-1};
    }
    for(bool mean:{false,true}) for(bool accumulate:{false,true}) {
      std::fill_n(a.data,n,.3f);std::fill_n(b.data,n,.3f);
      native.power(spectrum.data,mean ? grid.data : nullptr,.3f,a.data,n,accumulate);
      scalar.power(spectrum.data,mean ? grid.data : nullptr,.3f,b.data,n,accumulate);
      CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
      CHECK(native.score(a.data,weights.data,n)==scalar.score(b.data,weights.data,n));
      native.window(a.data,weights.data,n,.3f);scalar.window(b.data,weights.data,n,.3f);
      CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
      native.scale(a.data,weights.data,n,.5f,5);scalar.scale(b.data,weights.data,n,.5f,5);
      CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
      native.scale(a.data,nullptr,n,2,1);scalar.scale(b.data,nullptr,n,2,1);
      CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
    }
    // Adding one to 2^24 rounds back to 2^24. A separately reduced narrow tail
    // would change the answer once enough ones had been grouped together.
    std::fill_n(a.data,n,1.f);a.data[0]=16777216.f;std::fill_n(weights.data,n,1.f);
    CHECK(native.score(a.data,weights.data,n)==16777216.f);
    CHECK(native.score(a.data,weights.data,n)==scalar.score(a.data,weights.data,n));
    for(auto k:{native,scalar}) for(std::size_t bad=0;bad<n;++bad) {
      std::fill_n(a.data,n,1.f);a.data[bad]=INFINITY;
      rejects([&]{k.window(a.data,weights.data,n,0);});
      std::fill_n(a.data,n,1.f);a.data[bad]=NAN;
      rejects([&]{k.decode(a.data,SampleStorage::F32,b.data,n,0,0);});
      std::fill_n(a.data,n,1.f);a.data[bad]=1e30f;
      rejects([&]{k.scale(a.data,nullptr,n,1e20f,0);});
      const auto saved=spectrum.data[bad];spectrum.data[bad]={1e30f,0};
      rejects([&]{k.power(spectrum.data,nullptr,0,b.data,n,false);});spectrum.data[bad]=saved;
    }
  }
  std::cout<<"model kernels: exact scalar agreement, guard tails and non-finite checks passed\n";
}
int main() {try {
#ifdef NEO_FFT_TEST_HIGHWAY
  // Exercise only CPU-supported targets that were compiled into the core.
  // SetSupportedTargetsForTest also invalidates Highway's dispatch cache.
  const auto targets=hwy::SupportedAndGeneratedTargets();CHECK(!targets.empty());
  for(const auto target:targets) {
    hwy::SetSupportedTargetsForTest(target);
    CHECK(std::string(spatial_target(0))==hwy::TargetName(target));
    std::cout<<hwy::TargetName(target)<<": "<<std::flush;
    run_tests();
  }
  hwy::SetSupportedTargetsForTest(0);
#else
  run_tests();
#endif
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;} }
