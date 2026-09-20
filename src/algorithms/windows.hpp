#pragma once
#include "spectral/geometry.hpp"

namespace neo_fft {
struct AxisWindow {
  std::vector<float> analysis, synthesis;
};
AxisWindow fft3d_window(int block, int overlap, int type);
double dft_raw_window(int id, int j, int length, float beta);
struct DftWindow {
  std::vector<float> h;
  float wscale;
};
DftWindow dft_window(int block, int overlap, int mode, int spatial, int temporal, float sbeta, float tbeta);
} // namespace neo_fft
