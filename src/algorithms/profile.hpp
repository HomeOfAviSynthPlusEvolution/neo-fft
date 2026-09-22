#pragma once
#include "base/checked.hpp"
#include <array>

namespace neo_fft {
// Raw binary32 pairs; validation also applies to overridden axes.
struct DFTCurves {
  std::vector<float> shared, x, y, time;
  int system = 0;
  bool empty() const { return shared.empty() && x.empty() && y.empty() && time.empty(); }
};
void validate(const DFTCurves& curves);
// Returns calibrated primary values in [T,S,S/2+1] order.
std::vector<float> dft_profile(const DFTCurves& curves, int time, int size, float sigma, float divisor, int opt = 0);
} // namespace neo_fft
