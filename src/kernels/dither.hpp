#pragma once
#include "base/checked.hpp"
namespace neo_fft {
std::uint32_t dither_hash(std::uint32_t seed,std::uint32_t n,std::uint32_t p,std::uint32_t y,std::uint32_t x);
void dither_scalar(span2d::Plane<const float> source, span2d::Plane<std::uint8_t> dst,
                   int mode, int seed, int frame, int plane);
} // namespace neo_fft
