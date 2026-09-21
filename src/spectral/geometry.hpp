#pragma once
#include "base/checked.hpp"
#include <algorithm>

namespace neo_fft {
struct Axis {
  int length = 0, block = 0, overlap = 0, step = 0, count = 0, cover = 0, offset = 0;
};
inline Axis admit(Axis a) {
  require(a.offset >= 0 && a.cover >= std::int64_t(a.offset) + a.length && a.offset < a.length &&
              std::int64_t(a.cover) - a.offset - a.length < a.length,
          "selected-plane padding exceeds one-reflection domain");
  return a;
}
inline Axis fft3d_axis(int length, int block, int overlap) {
  dimension(length);
  require(block >= 2 && length >= block, "FFT3D block must be >=2 and fit selected plane");
  if (overlap < 0)
    overlap = block / 3;
  require(overlap <= block / 2, "FFT3D overlap exceeds half block");
  const int step = block - overlap;
  const int count = dimension(ceil_div(std::int64_t(length) - overlap, step) + 2);
  const int cover = dimension(std::int64_t(count) * step + overlap);
  return admit({length, block, overlap, step, count, cover, step});
}
inline Axis dft_axis(int length, int block, int mode, int overlap) {
  dimension(length);
  dimension(block);
  require(mode == 0 || mode == 1, "DFTTest smode must be 0 or 1");
  if (mode == 0) {
    require(block % 2 == 1, "DFTTest center mode requires odd sbsize");
    const int c = block / 2;
    return admit({length, block, 0, 1, length, dimension(std::int64_t(length) + 2LL * c), c});
  }
  require(overlap >= 0 && overlap < block, "DFTTest sosize outside block");
  const int step = block - overlap;
  require(overlap <= block / 2 || block % step == 0, "DFTTest heavy overlap requires divisible step");
  const int cover = dimension(std::int64_t(block) * ceil_div(length, block) + 2LL * std::max(step, overlap));
  Axis a = admit({length, block, overlap, step, (cover - block) / step + 1, cover, (cover - length) / 2});
  // Every spatial phase in the crop has all of its contributing windows.
  require(a.offset >= block - step && std::int64_t(a.offset) + length <= std::int64_t(a.count) * step,
          "DFTTest visible crop lacks complete phase coverage");
  return a;
}
inline int reflect(int cover_coordinate, const Axis& a) {
  const std::int64_t j = std::int64_t(cover_coordinate) - a.offset;
  return static_cast<int>(j < 0 ? -j : j < a.length ? j : 2LL * a.length - 2 - j);
}
struct Geometry {
  Axis x, y;
  Geometry(Axis x_axis, Axis y_axis) : x(x_axis), y(y_axis) {
    plane_extent<float>(x.cover, y.cover, std::ptrdiff_t(mul_size(x.cover, sizeof(float))));
    mul_size(std::size_t(x.count), std::size_t(y.count));
  }
};
} // namespace neo_fft
