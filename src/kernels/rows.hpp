#pragma once
#include "base/checked.hpp"
namespace neo_fft {
using CopyRow=void (*)(const void*,void*,std::size_t);
using DitherNoise=void (*)(float*,std::size_t,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t);
void copy_row_scalar(const void* source,void* dst,std::size_t bytes);
void dither_noise_scalar(float* dst,std::size_t count,std::uint32_t seed,std::uint32_t frame,std::uint32_t plane,std::uint32_t y,std::uint32_t x);
CopyRow select_copy_row(int opt);
DitherNoise select_dither_noise(int opt);
} // namespace neo_fft
