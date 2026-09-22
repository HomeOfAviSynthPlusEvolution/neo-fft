#pragma once
#include "kernels/spectral.hpp"
#include <array>

namespace neo_fft {
struct EnhancementConfig {
  float sharpen = 0, scutoff = .3f, svr = 1, smin = 4, smax = 20;
  float dehalo = 0, hr = 2, ht = 50;
};
void validate(const EnhancementConfig& c);
std::vector<float> fft3d_profile(int width, int height, const std::array<float, 4>& sigmas, int opt = 0);
struct EnhancementTables {
  std::vector<float> sharpen, halo;
  float a = 0, b = 0, c = 0;
  Enhancement view(const EnhancementConfig& config) const;
};
EnhancementTables enhancement_tables(int width, int height, float format_scale, const EnhancementConfig& c, int opt = 0);
} // namespace neo_fft
