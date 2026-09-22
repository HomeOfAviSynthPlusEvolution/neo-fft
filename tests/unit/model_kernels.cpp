#include "kernels/model.hpp"
#include "base/checked.hpp"
#include "../test.hpp"
#include <cstring>
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
  for(std::size_t i=0;i<n;++i) src.data[i]=T(i%251);
  for(float base:{0.f,128.f}) for(float scale:{1.f,255.f,1.f/256}) {
    native.decode(src.data,sample_storage<T>,a.data,n,base,scale);
    model_scalar().decode(src.data,sample_storage<T>,b.data,n,base,scale);
    CHECK(std::memcmp(a.data,b.data,n*sizeof(float))==0);
  }
}
int main() {try {
  const auto native=select_model(0);const auto scalar=model_scalar();
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
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;} }
