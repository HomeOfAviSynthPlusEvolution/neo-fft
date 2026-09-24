#pragma once
#include <cstddef>
#include <cstdint>

namespace neo_fft {

// Inspect exactly count contiguous samples; zero count permits nullptr.
// Throws the same error as finite(), without changing any samples.
using ValidateFiniteFn = void (*)(const float* src, std::size_t count);
using GatherFft3dFn = void (*)(const float* src, const float* wx_a, float wy, float* blk, int count) noexcept;
// Gather one square spatial block from each temporal slice. Strides are in floats;
// window and output contain size*size*temporal tightly packed samples. The caller
// supplies positive dimensions and valid, non-overlapping input/output storage.
using GatherDfttestFn = void (*)(const float* src, std::ptrdiff_t row_stride, std::size_t slice_stride,
                                const float* window, float* block, int size, int temporal) noexcept;
using ScatterFft3dBlockFn = void (*)(const float* inv, const float* wx_s, float* r_row, int count) noexcept;
using ScatterFft3dRowFn = void (*)(const float* r_ptr, float wy, float* a_ptr, int count) noexcept;
using ScatterDfttestFn = void (*)(const float* inv, const float* h_syn, float* acc_row, int count) noexcept;
using StoreOutputFloatFn = void (*)(const float* a_ptr, float* dst_row, int count, bool fft3d, float scale) noexcept;
using StoreOutputUint8Fn = void (*)(const float* a_ptr, std::uint8_t* dst_row, int count, bool fft3d, float base,
                                    float scale, float peak) noexcept;
using StoreOutputUint16Fn = void (*)(const float* a_ptr, std::uint16_t* dst_row, int count, bool fft3d, float base,
                                     float scale, float peak) noexcept;

struct SpatialKernels {
  ValidateFiniteFn validate_finite = nullptr;
  GatherFft3dFn gather_fft3d = nullptr;
  GatherDfttestFn gather_dfttest = nullptr;
  ScatterFft3dBlockFn scatter_fft3d_block = nullptr;
  ScatterFft3dRowFn scatter_fft3d_row = nullptr;
  ScatterDfttestFn scatter_dfttest = nullptr;
  StoreOutputFloatFn store_output_float = nullptr;
  StoreOutputUint8Fn store_output_uint8 = nullptr;
  StoreOutputUint16Fn store_output_uint16 = nullptr;
};

SpatialKernels select_spatial(int opt);
const SpatialKernels& spatial_scalar() noexcept;
const char* spatial_target(int opt);

} // namespace neo_fft
