#include "spectral/fft_backend.hpp"
#include "base/checked.hpp"
#include "spectral/codelets/fft16x16.hpp"
#include "spectral/codelets/fft8x8.hpp"
#include "spectral/codelets/fft32x32.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
class TlsScratchPool {
public:
  static constexpr std::size_t CAPACITY = 1024 * 1024; // 1 MB

  TlsScratchPool() = default;
  ~TlsScratchPool() {
    if (buffer_) {
#if defined(_WIN32)
      _aligned_free(buffer_);
#else
      std::free(buffer_);
#endif
      buffer_ = nullptr;
    }
  }

  TlsScratchPool(const TlsScratchPool&) = delete;
  TlsScratchPool& operator=(const TlsScratchPool&) = delete;

  void* allocate(std::size_t size) {
    if (active_) {
      if (!buffer_) {
#if defined(_WIN32)
        buffer_ = static_cast<char*>(_aligned_malloc(CAPACITY, 64));
#else
        void* ptr = nullptr;
        if (posix_memalign(&ptr, 64, CAPACITY) != 0) {
          ptr = std::malloc(CAPACITY);
        }
        buffer_ = static_cast<char*>(ptr);
#endif
      }
      if (buffer_) {
        const std::size_t aligned_size = (size + 63) & ~std::size_t(63);
        if (offset_ + aligned_size <= CAPACITY) {
          char* ptr = buffer_ + offset_;
          offset_ += aligned_size;
          ++active_allocs_;
          last_ptr_ = ptr;
          last_size_ = aligned_size;
          return ptr;
        }
      }
    }
    return std::malloc(size);
  }

  void deallocate(void* ptr) noexcept {
    if (!ptr)
      return;
    if (buffer_ && ptr >= buffer_ && ptr < buffer_ + CAPACITY) {
      if (active_allocs_ > 0) {
        --active_allocs_;
      }
      if (active_allocs_ == 0) {
        offset_ = 0;
        last_ptr_ = nullptr;
        last_size_ = 0;
      } else if (ptr == last_ptr_) {
        offset_ -= last_size_;
        last_ptr_ = nullptr;
        last_size_ = 0;
      }
      return;
    }
    std::free(ptr);
  }

  bool is_active() const noexcept {
    return active_;
  }

  void set_active(bool active) noexcept {
    active_ = active;
  }

  std::size_t active_allocs() const noexcept {
    return active_allocs_;
  }

  void reset() noexcept {
    if (active_allocs_ == 0) {
      offset_ = 0;
      last_ptr_ = nullptr;
      last_size_ = 0;
    }
  }

private:
  char* buffer_ = nullptr;
  std::size_t offset_ = 0;
  std::size_t active_allocs_ = 0;
  char* last_ptr_ = nullptr;
  std::size_t last_size_ = 0;
  bool active_ = false;
};

thread_local TlsScratchPool g_scratch_pool;

inline void* pocketfft_scratch_malloc(std::size_t size) {
  return g_scratch_pool.allocate(size);
}

inline void pocketfft_scratch_free(void* ptr) noexcept {
  g_scratch_pool.deallocate(ptr);
}

class ScratchScope {
public:
  ScratchScope() {
    g_scratch_pool.set_active(true);
  }
  ~ScratchScope() {
    g_scratch_pool.set_active(false);
    g_scratch_pool.reset();
  }
  ScratchScope(const ScratchScope&) = delete;
  ScratchScope& operator=(const ScratchScope&) = delete;
};

class SuspendScratchScope {
public:
  SuspendScratchScope() noexcept : prev_(g_scratch_pool.is_active()) {
    g_scratch_pool.set_active(false);
  }
  ~SuspendScratchScope() noexcept {
    g_scratch_pool.set_active(prev_);
  }
  SuspendScratchScope(const SuspendScratchScope&) = delete;
  SuspendScratchScope& operator=(const SuspendScratchScope&) = delete;

private:
  bool prev_;
};
} // namespace

// Built multiple times with separate dependency namespaces. No shared
// cache or worker pool exists between targets.
#if defined(NEO_FFT_FFT_SCALAR)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_scalar_c90e55b3
#ifndef POCKETFFT_NO_VECTORS
#define POCKETFFT_NO_VECTORS
#endif
#define BACKEND_FN scalar_fft
#define BACKEND_NAME "pocketfft-scalar"
#elif defined(NEO_FFT_FFT_SSE2)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_sse2_c90e55b3
#define BACKEND_FN sse2_fft
#define BACKEND_NAME "pocketfft-sse2"
#elif defined(NEO_FFT_FFT_AVX2)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_avx2_c90e55b3
#define BACKEND_FN avx2_fft
#define BACKEND_NAME "pocketfft-avx2"
#elif defined(NEO_FFT_FFT_AVX512)
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_avx512_c90e55b3
#define BACKEND_FN avx512_fft
#define BACKEND_NAME "pocketfft-avx512"
#else
#define POCKETFFT_NAMESPACE neo_fft_pocketfft_native_c90e55b3
#define BACKEND_FN native_target_fft
#define BACKEND_NAME "pocketfft-native"
#endif
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_CACHE_SIZE 16

#define malloc(sz) pocketfft_scratch_malloc(sz)
#define free(ptr) pocketfft_scratch_free(ptr)
#include <pocketfft_hdronly.h>
#undef malloc
#undef free

namespace POCKETFFT_NAMESPACE::detail {

template <typename T>
inline std::shared_ptr<T> get_plan_persistent(size_t length) {
  static std::mutex mut;
  static std::array<std::pair<size_t,std::shared_ptr<T>>,POCKETFFT_CACHE_SIZE> cache{};
  static size_t next=0;
  // Weak TLS entries skip the shared-cache lock for a thread's active axes
  // without extending the lifetime of plans evicted from the bounded cache.
  static thread_local std::array<std::pair<size_t,std::weak_ptr<T>>,4> recent{};
  static thread_local size_t next_recent=0;
  for(const auto& entry:recent)
    if(entry.first==length) if(auto plan=entry.second.lock()) return plan;
  const auto remember=[&](std::shared_ptr<T> plan) {
    for(auto& entry:recent) if(entry.first==length) {entry.second=plan;return plan;}
    recent[next_recent]={length,plan};
    next_recent=(next_recent+1)%recent.size();
    return plan;
  };
  {
    std::lock_guard<std::mutex> lock(mut);
    for(const auto& entry:cache) if(entry.first==length && entry.second) return remember(entry.second);
  }
  SuspendScratchScope suspend;
  auto plan = std::make_shared<T>(length);
  {
    std::lock_guard<std::mutex> lock(mut);
    for(const auto& entry:cache) if(entry.first==length && entry.second) return remember(entry.second);
    cache[next]={length,plan};
    next=(next+1)%cache.size();
  }
  return remember(std::move(plan));
}

template <>
inline std::shared_ptr<pocketfft_r<float>> get_plan<pocketfft_r<float>>(size_t length) {
  return get_plan_persistent<pocketfft_r<float>>(length);
}

template <>
inline std::shared_ptr<pocketfft_c<float>> get_plan<pocketfft_c<float>>(size_t length) {
  return get_plan_persistent<pocketfft_c<float>>(length);
}

template <>
inline std::shared_ptr<pocketfft_r<double>> get_plan<pocketfft_r<double>>(size_t length) {
  return get_plan_persistent<pocketfft_r<double>>(length);
}

template <>
inline std::shared_ptr<pocketfft_c<double>> get_plan<pocketfft_c<double>>(size_t length) {
  return get_plan_persistent<pocketfft_c<double>>(length);
}

} // namespace POCKETFFT_NAMESPACE::detail

namespace neo_fft::detail {
namespace {
namespace pf = POCKETFFT_NAMESPACE;

void r2c(int height, int width, const float* in, std::size_t in_row_stride, std::complex<float>* out,
         std::size_t out_row_stride) {
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::fft32x32_r2c(in, in_row_stride, out, out_row_stride);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::fft16x16_r2c(in, in_row_stride, out, out_row_stride);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::fft8x8_r2c(in, in_row_stride, out, out_row_stride);
    return;
  }
#endif
  ScratchScope scope;
  const pf::shape_t shape{std::size_t(height), std::size_t(width)}, axes{0, 1};
  const pf::stride_t rs{std::ptrdiff_t(in_row_stride * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void c2r(int height, int width, const std::complex<float>* in, std::size_t in_row_stride, float* out,
         std::size_t out_row_stride, float fct) {
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::fft32x32_c2r(in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::fft16x16_c2r(in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::fft8x8_c2r(in, in_row_stride, out, out_row_stride, fct);
    return;
  }
#endif
  ScratchScope scope;
  const pf::shape_t shape{std::size_t(height), std::size_t(width)}, axes{0, 1};
  const pf::stride_t rs{std::ptrdiff_t(out_row_stride * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}

void batch_r2c(std::size_t batch, int height, int width, const float* in, std::size_t in_dist,
               std::size_t in_row_stride, std::complex<float>* out, std::size_t out_dist,
               std::size_t out_row_stride) {
  if (batch == 0)
    return;
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::batch_fft32x32_r2c(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::batch_fft16x16_r2c(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::batch_fft8x8_r2c(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride);
    return;
  }
#endif
  if (batch == 1) {
    r2c(height, width, in, in_row_stride, out, out_row_stride);
    return;
  }
  ScratchScope scope;
  const pf::shape_t shape{batch, std::size_t(height), std::size_t(width)}, axes{1, 2};
  const pf::stride_t rs{std::ptrdiff_t(in_dist * sizeof(float)), std::ptrdiff_t(in_row_stride * sizeof(float)),
                        sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(out_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void batch_c2r(std::size_t batch, int height, int width, const std::complex<float>* in, std::size_t in_dist,
               std::size_t in_row_stride, float* out, std::size_t out_dist, std::size_t out_row_stride, float fct) {
  if (batch == 0)
    return;
#if defined(NEO_FFT_HAS_AVX2_CODELET)
  if (height == 32 && width == 32) {
    codelet::batch_fft32x32_c2r(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride, fct);
    return;
  }
  if (height == 16 && width == 16) {
    codelet::batch_fft16x16_c2r(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride, fct);
    return;
  }
  if (height == 8 && width == 8) {
    codelet::batch_fft8x8_c2r(batch, in, in_dist, in_row_stride, out, out_dist, out_row_stride, fct);
    return;
  }
#endif
  if (batch == 1) {
    c2r(height, width, in, in_row_stride, out, out_row_stride, fct);
    return;
  }
  ScratchScope scope;
  const pf::shape_t shape{batch, std::size_t(height), std::size_t(width)}, axes{1, 2};
  const pf::stride_t rs{std::ptrdiff_t(out_dist * sizeof(float)), std::ptrdiff_t(out_row_stride * sizeof(float)),
                        sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_dist * sizeof(std::complex<float>)),
                        std::ptrdiff_t(in_row_stride * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}

inline std::ptrdiff_t checked_stride(int rows,int columns,std::size_t bytes) {
  const auto result=mul_size(mul_size(std::size_t(rows),std::size_t(columns)),bytes);
  require(result<=std::size_t(PTRDIFF_MAX),"FFT stride exceeds ptrdiff_t");
  return static_cast<std::ptrdiff_t>(result);
}
void r2c_3d(int depth, int height, int width, const float* in, std::complex<float>* out) {
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{0, 1, 2};
  const pf::stride_t rs{checked_stride(height,width,sizeof(float)), std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{checked_stride(height,k,sizeof(std::complex<float>)), std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void c2r_3d(int depth, int height, int width, const std::complex<float>* in, float* out, float fct) {
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{0, 1, 2};
  const pf::stride_t rs{checked_stride(height,width,sizeof(float)), std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{checked_stride(height,k,sizeof(std::complex<float>)), std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}

// Dense small-volume C2C axes avoid the general multidimensional iterator.
// Keep PocketFFT's butterflies and axis order: width real FFT, time, height.
#ifndef POCKETFFT_NO_VECTORS
template<int N,int Inner,int Outer>
void dense_axis(std::size_t batch,const std::complex<float>* in,std::size_t in_dist,
                std::complex<float>* out,std::size_t out_dist,bool forward) {
  namespace pd=pf::detail;
  constexpr auto lanes=pd::VLEN<float>::val;
  using Vec=pd::vtype_t<float>;
  const auto plan=pd::get_plan<pd::pocketfft_c<float>>(N);
  const auto lines=batch*Outer*Inner;
  std::size_t line=0;
  for(;lines-line>=lanes;line+=lanes) {
    std::array<std::size_t,lanes> src{},dst{};
    for(std::size_t lane=0;lane<lanes;++lane) {
      const auto i=line+lane,b=i/(Outer*Inner),within=i%(Outer*Inner);
      const auto offset=(within/Inner)*N*Inner+within%Inner;
      src[lane]=b*in_dist+offset;dst[lane]=b*out_dist+offset;
    }
    pd::cmplx<Vec> values[N];
    for(int n=0;n<N;++n)for(std::size_t lane=0;lane<lanes;++lane) {
      const auto v=in[src[lane]+n*Inner];
      values[n].r[lane]=v.real();values[n].i[lane]=v.imag();
    }
    plan->exec(values,1.0f,forward);
    for(int n=0;n<N;++n)for(std::size_t lane=0;lane<lanes;++lane)
      out[dst[lane]+n*Inner]={values[n].r[lane],values[n].i[lane]};
  }
  for(;line<lines;++line) {
    const auto b=line/(Outer*Inner),within=line%(Outer*Inner);
    const auto offset=(within/Inner)*N*Inner+within%Inner;
    pd::cmplx<float> values[N];
    for(int n=0;n<N;++n) {
      const auto v=in[b*in_dist+offset+n*Inner];values[n].Set(v.real(),v.imag());
    }
    plan->exec(values,1.0f,forward);
    for(int n=0;n<N;++n)out[b*out_dist+offset+n*Inner]={values[n].r,values[n].i};
  }
}

template<int T,int S>
void dense_axes(std::size_t batch,const std::complex<float>* in,std::size_t in_dist,
                std::complex<float>* out,std::size_t out_dist,bool forward,bool spatial) {
  constexpr int K=S/2+1;
  dense_axis<T,S*K,1>(batch,in,in_dist,out,out_dist,forward);
  if(spatial)dense_axis<S,K,T>(batch,out,out_dist,out,out_dist,forward);
}

bool dense_shape(int depth,int height,int width) {
  return (depth==3 || depth==5) && (height==12 || height==16) && width==height;
}
void dense_axes_dispatch(std::size_t batch,int depth,int size,const std::complex<float>* in,
                         std::size_t in_dist,std::complex<float>* out,std::size_t out_dist,bool forward,bool spatial=true) {
  if(size==12) {
    if(depth==3)dense_axes<3,12>(batch,in,in_dist,out,out_dist,forward,spatial);
    else dense_axes<5,12>(batch,in,in_dist,out,out_dist,forward,spatial);
  } else {
    if(depth==3)dense_axes<3,16>(batch,in,in_dist,out,out_dist,forward,spatial);
    else dense_axes<5,16>(batch,in,in_dist,out,out_dist,forward,spatial);
  }
}
#endif

void batch_r2c_3d(std::size_t batch, int depth, int height, int width, const float* in, std::size_t in_dist,
                  std::complex<float>* out, std::size_t out_dist) {
  if (batch == 0)
    return;
  if (batch == 1) {
    r2c_3d(depth, height, width, in, out);
    return;
  }
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{batch, std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{1, 2, 3};
  const pf::stride_t rs{std::ptrdiff_t(in_dist * sizeof(float)), checked_stride(height,width,sizeof(float)),
                        std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(out_dist * sizeof(std::complex<float>)),
                        checked_stride(height,k,sizeof(std::complex<float>)),
                        std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
#ifndef POCKETFFT_NO_VECTORS
  if(dense_shape(depth,height,width)) {
#if defined(NEO_FFT_HAS_AVX2_CODELET)
    if(width==16) {
      // Spatial codelets first, then the temporal axis. This is used only by
      // grouped callers; single-volume threshold-sensitive paths stay generic.
      for(std::size_t b=0;b<batch;++b)
        codelet::batch_fft16x16_r2c(depth,in+b*in_dist,256,16,out+b*out_dist,144,9);
      dense_axes_dispatch(batch,depth,height,out,out_dist,out,out_dist,true,false);
      return;
    }
#endif
    pf::r2c(shape,rs,ss,std::size_t(3),true,in,out,1.0f,1);
    dense_axes_dispatch(batch,depth,height,out,out_dist,out,out_dist,true);
    return;
  }
#endif
  pf::r2c(shape, rs, ss, axes, true, in, out, 1.0f, 1);
}

void batch_c2r_3d(std::size_t batch, int depth, int height, int width, const std::complex<float>* in, std::size_t in_dist,
                  float* out, std::size_t out_dist, float fct) {
  if (batch == 0)
    return;
  if (batch == 1) {
    c2r_3d(depth, height, width, in, out, fct);
    return;
  }
  ScratchScope scope;
  const int k = width / 2 + 1;
  const pf::shape_t shape{batch, std::size_t(depth), std::size_t(height), std::size_t(width)}, axes{1, 2, 3};
  const pf::stride_t rs{std::ptrdiff_t(out_dist * sizeof(float)), checked_stride(height,width,sizeof(float)),
                        std::ptrdiff_t(width * sizeof(float)), sizeof(float)};
  const pf::stride_t ss{std::ptrdiff_t(in_dist * sizeof(std::complex<float>)),
                        checked_stride(height,k,sizeof(std::complex<float>)),
                        std::ptrdiff_t(k * sizeof(std::complex<float>)), sizeof(std::complex<float>)};
#ifndef POCKETFFT_NO_VECTORS
  if(dense_shape(depth,height,width)) {
    const auto bins=std::size_t(depth)*height*k;
    pf::detail::arr<std::complex<float>> temporary(mul_size(batch,bins));
#if defined(NEO_FFT_HAS_AVX2_CODELET)
    if(width==16) {
      dense_axes_dispatch(batch,depth,height,in,in_dist,temporary.data(),bins,false,false);
      for(std::size_t b=0;b<batch;++b)
        codelet::batch_fft16x16_c2r(depth,temporary.data()+b*bins,144,9,out+b*out_dist,256,16,fct);
      return;
    }
#endif
    dense_axes_dispatch(batch,depth,height,in,in_dist,temporary.data(),bins,false);
    auto dense_stride=ss;dense_stride[0]=std::ptrdiff_t(bins*sizeof(std::complex<float>));
    pf::c2r(shape,dense_stride,rs,std::size_t(3),false,temporary.data(),out,fct,1);
    return;
  }
#endif
  pf::c2r(shape, ss, rs, axes, false, in, out, fct, 1);
}
} // namespace

const FftBackend& BACKEND_FN() noexcept {
  static const FftBackend backend{int(pf::detail::VLEN<float>::val), BACKEND_NAME, r2c, c2r, batch_r2c, batch_c2r,
                                  r2c_3d, c2r_3d, batch_r2c_3d, batch_c2r_3d};
  return backend;
}

} // namespace neo_fft::detail
