#include "kernels/spatial.hpp"
#include "base/checked.hpp"
#include <algorithm>

namespace neo_fft {
namespace {

void gather_fft3d_scalar(const float* src, const float* wx_a, float wy, float* blk, int count) noexcept {
  for (int x = 0; x < count; ++x) {
    blk[x] = (src[x] * wy) * wx_a[x];
  }
}

void gather_dfttest_scalar(const float* src, const float* h_row, float* blk, int count) noexcept {
  for (int x = 0; x < count; ++x) {
    blk[x] = src[x] * h_row[x];
  }
}

void scatter_fft3d_block_scalar(const float* inv, const float* wx_s, float* r_row, int count) noexcept {
  for (int x = 0; x < count; ++x) {
    r_row[x] += inv[x] * wx_s[x];
  }
}

void scatter_fft3d_row_scalar(const float* r_ptr, float wy, float* a_ptr, int count) noexcept {
  for (int x = 0; x < count; ++x) {
    a_ptr[x] += r_ptr[x] * wy;
  }
}

void scatter_dfttest_scalar(const float* inv, const float* h_syn, float* acc_row, int count) noexcept {
  for (int x = 0; x < count; ++x) {
    acc_row[x] += inv[x] * h_syn[x];
  }
}

void store_output_float_scalar(const float* a_ptr, float* dst_row, int count, bool fft3d, float scale) noexcept {
  if (fft3d) {
    for (int x = 0; x < count; ++x) {
      dst_row[x] = std::clamp(a_ptr[x], 0.0f, 1.0f);
    }
  } else {
    for (int x = 0; x < count; ++x) {
      dst_row[x] = a_ptr[x] * scale;
    }
  }
}

void store_output_uint8_scalar(const float* a_ptr, std::uint8_t* dst_row, int count, bool fft3d, float base,
                               float scale, float peak) noexcept {
  if (fft3d) {
    for (int x = 0; x < count; ++x) {
      const float v = (a_ptr[x] + 0.5f) + base;
      dst_row[x] = static_cast<std::uint8_t>(std::clamp(v, 0.0f, peak));
    }
  } else {
    for (int x = 0; x < count; ++x) {
      const float v = (a_ptr[x] * scale) + 0.5f;
      dst_row[x] = static_cast<std::uint8_t>(std::clamp(v, 0.0f, peak));
    }
  }
}

void store_output_uint16_scalar(const float* a_ptr, std::uint16_t* dst_row, int count, bool fft3d, float base,
                                float scale, float peak) noexcept {
  if (fft3d) {
    for (int x = 0; x < count; ++x) {
      const float v = (a_ptr[x] + 0.5f) + base;
      dst_row[x] = static_cast<std::uint16_t>(std::clamp(v, 0.0f, peak));
    }
  } else {
    for (int x = 0; x < count; ++x) {
      const float v = (a_ptr[x] * scale) + 0.5f;
      dst_row[x] = static_cast<std::uint16_t>(std::clamp(v, 0.0f, peak));
    }
  }
}

} // namespace

const SpatialKernels& spatial_scalar() noexcept {
  static const SpatialKernels kScalar{
      gather_fft3d_scalar,
      gather_dfttest_scalar,
      scatter_fft3d_block_scalar,
      scatter_fft3d_row_scalar,
      scatter_dfttest_scalar,
      store_output_float_scalar,
      store_output_uint8_scalar,
      store_output_uint16_scalar,
  };
  return kScalar;
}

#if !NEO_FFT_ENABLE_HIGHWAY
SpatialKernels select_spatial(int opt) {
  return spatial_scalar();
}
const char* spatial_target(int opt) {
  return "scalar (Highway disabled)";
}
#endif

} // namespace neo_fft
