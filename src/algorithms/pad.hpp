#pragma once
#include "algorithms/plan.hpp"
#include <dualsynth/span2d.hpp>
#include <algorithm>

namespace neo_fft {

inline int reflect(int j, int L) noexcept {
  if (j < 0) return static_cast<int>(-std::int64_t(j));
  if (j >= L) return static_cast<int>(2LL * L - 2 - j);
  return j;
}

template <class T>
void pad_source(span2d::Plane<const T> src, span2d::Plane<float> padded, const Geometry& geom, SampleFormat format, Algorithm alg) {
  float base = 0.0f;
  float scale = 1.0f;

  if (alg == Algorithm::FFT3D) {
    base = (!format.floating && format.chroma) ? float(1 << (format.bits - 1)) : 0.0f;
  } else {
    scale = format.floating ? 255.0f : 1.0f / float(1 << (format.bits - 8));
  }

  const int L_w = src.width();
  const int L_h = src.height();
  const int P_w = geom.x.cover;
  const int P_h = geom.y.cover;
  const int d_x = geom.x.offset;
  const int d_y = geom.y.offset;

  for (int cy = 0; cy < P_h; ++cy) {
    int y = cy - d_y;
    int src_y = reflect(y, L_h);
    const T* src_row = src.row_ptr(src_y);
    float* pad_row = padded.row_ptr(cy);
    for (int cx = 0; cx < P_w; ++cx) {
      int x = cx - d_x;
      int src_x = reflect(x, L_w);
      float val = float(src_row[src_x]);
      if (alg == Algorithm::FFT3D) {
        val = val - base;
      } else {
        val = val * scale;
      }
      pad_row[cx] = val;
    }
  }
}

} // namespace neo_fft
