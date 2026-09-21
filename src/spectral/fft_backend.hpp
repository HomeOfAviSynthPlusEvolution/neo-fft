#pragma once
#include <complex>
#include <cstddef>

namespace neo_fft {
enum class FftProfile { scalar, native, sse2, avx2, avx512 };

struct FftBackend {
  int lanes;
  const char* name;
  void (*r2c)(int height, int width, const float* in, std::size_t in_row_stride, std::complex<float>* out,
              std::size_t out_row_stride);
  void (*c2r)(int height, int width, const std::complex<float>* in, std::size_t in_row_stride, float* out,
              std::size_t out_row_stride, float fct);
  void (*batch_r2c)(std::size_t batch, int height, int width, const float* in, std::size_t in_dist,
                    std::size_t in_row_stride, std::complex<float>* out, std::size_t out_dist,
                    std::size_t out_row_stride);
  void (*batch_c2r)(std::size_t batch, int height, int width, const std::complex<float>* in, std::size_t in_dist,
                    std::size_t in_row_stride, float* out, std::size_t out_dist,
                    std::size_t out_row_stride, float fct);
};

namespace detail {
const FftBackend& scalar_fft() noexcept;
#if NEO_FFT_FFT_X86_TARGETS
const FftBackend& sse2_fft() noexcept;
const FftBackend& avx2_fft() noexcept;
const FftBackend& avx512_fft() noexcept;
#else
const FftBackend& native_target_fft() noexcept;
#endif
const FftBackend& native_fft() noexcept;
} // namespace detail

int fft_lanes(FftProfile profile = FftProfile::native) noexcept;
const char* fft_profile_name(FftProfile profile = FftProfile::native) noexcept;
const char* fft_backend_name(FftProfile profile = FftProfile::native) noexcept;
const FftBackend& backend_by_profile(FftProfile profile) noexcept;
bool fft_profile_supported(FftProfile profile = FftProfile::native) noexcept;

} // namespace neo_fft
