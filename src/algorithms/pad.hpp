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

// Admitted geometry contains the full source and uses one reflection.
// Decode and validate each source sample once, then copy the reflected borders.
template <class T>
void pad_decoded_source(span2d::Plane<const T> src, span2d::Plane<float> padded,
                        const Geometry& geom, float base, float scale, const ModelKernels& model) {
  const int width = src.width(), height = src.height();
  const int dx = geom.x.offset, dy = geom.y.offset;
  for (int y = 0; y < height; ++y) {
    float* row = padded.row_ptr(dy + y);
    model.decode(src.row_ptr(y), sample_storage<T>, row + dx, std::size_t(width), base, scale);
    for (int x = 0; x < dx; ++x) row[x] = row[dx + reflect(x - dx, width)];
    for (int x = dx + width; x < geom.x.cover; ++x) row[x] = row[dx + reflect(x - dx, width)];
  }
  for (int y = 0; y < dy; ++y)
    std::copy_n(padded.row_ptr(dy + reflect(y - dy, height)), geom.x.cover, padded.row_ptr(y));
  for (int y = dy + height; y < geom.y.cover; ++y)
    std::copy_n(padded.row_ptr(dy + reflect(y - dy, height)), geom.x.cover, padded.row_ptr(y));
}

template <class T>
void pad_fft3d_source(span2d::Plane<const T> src, span2d::Plane<float> padded,
                      const Geometry& geom, SampleFormat format, const ModelKernels& model) {
  const float base = !format.floating && format.chroma ? float(1 << (format.bits - 1)) : 0.0f;
  pad_decoded_source(src,padded,geom,base,1.0f,model);
}

template <class T>
void pad_dfttest_source(span2d::Plane<const T> src, span2d::Plane<float> padded,
                       const Geometry& geom, SampleFormat format, const ModelKernels& model) {
  const float scale=format.floating ? 255.0f : 1.0f/float(1 << (format.bits-8));
  pad_decoded_source(src,padded,geom,0.0f,scale,model);
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
