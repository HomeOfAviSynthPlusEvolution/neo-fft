#pragma once
#include "test.hpp"
#include <vector>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
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
