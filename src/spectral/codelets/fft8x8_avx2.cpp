#include "spectral/codelets/fft8x8.hpp"

#if defined(__AVX2__)
#include <immintrin.h>
#include <cmath>

namespace neo_fft::codelet {
namespace {

constexpr float SQRT2_2 = 0.70710678118654752440f;

inline void transpose8x8_ps(__m256& row0, __m256& row1, __m256& row2, __m256& row3,
                            __m256& row4, __m256& row5, __m256& row6, __m256& row7) noexcept {
  __m256 t0 = _mm256_unpacklo_ps(row0, row1);
  __m256 t1 = _mm256_unpackhi_ps(row0, row1);
  __m256 t2 = _mm256_unpacklo_ps(row2, row3);
  __m256 t3 = _mm256_unpackhi_ps(row2, row3);
  __m256 t4 = _mm256_unpacklo_ps(row4, row5);
  __m256 t5 = _mm256_unpackhi_ps(row4, row5);
  __m256 t6 = _mm256_unpacklo_ps(row6, row7);
  __m256 t7 = _mm256_unpackhi_ps(row6, row7);

  __m256 sh0 = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1, 0, 1, 0));
  __m256 sh1 = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 2, 3, 2));
  __m256 sh2 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));
  __m256 sh3 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));
  __m256 sh4 = _mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(1, 0, 1, 0));
  __m256 sh5 = _mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(3, 2, 3, 2));
  __m256 sh6 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0));
  __m256 sh7 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2));

  row0 = _mm256_permute2f128_ps(sh0, sh4, 0x20);
  row1 = _mm256_permute2f128_ps(sh1, sh5, 0x20);
  row2 = _mm256_permute2f128_ps(sh2, sh6, 0x20);
  row3 = _mm256_permute2f128_ps(sh3, sh7, 0x20);
  row4 = _mm256_permute2f128_ps(sh0, sh4, 0x31);
  row5 = _mm256_permute2f128_ps(sh1, sh5, 0x31);
  row6 = _mm256_permute2f128_ps(sh2, sh6, 0x31);
  row7 = _mm256_permute2f128_ps(sh3, sh7, 0x31);
}

// 8-way 1D 8-point complex FFT
template <bool FORWARD>
inline void fft8_8way_avx2(__m256* r, __m256* i) noexcept {
  const __m256 v_sqrt2_2 = _mm256_set1_ps(SQRT2_2);
  const __m256 v_zero = _mm256_setzero_ps();

  // Bit reversal permutation: 0, 4, 2, 6, 1, 5, 3, 7
  __m256 a0_r = r[0], a0_i = i[0];
  __m256 a1_r = r[4], a1_i = i[4];
  __m256 a2_r = r[2], a2_i = i[2];
  __m256 a3_r = r[6], a3_i = i[6];
  __m256 a4_r = r[1], a4_i = i[1];
  __m256 a5_r = r[5], a5_i = i[5];
  __m256 a6_r = r[3], a6_i = i[3];
  __m256 a7_r = r[7], a7_i = i[7];

  #define B2_8(x_r, x_i, y_r, y_i) { \
    __m256 ur = x_r, ui = x_i; \
    __m256 vr = y_r, vi = y_i; \
    x_r = _mm256_add_ps(ur, vr); \
    x_i = _mm256_add_ps(ui, vi); \
    y_r = _mm256_sub_ps(ur, vr); \
    y_i = _mm256_sub_ps(ui, vi); \
  }

  // Stage 1: Length 2
  B2_8(a0_r, a0_i, a1_r, a1_i);
  B2_8(a2_r, a2_i, a3_r, a3_i);
  B2_8(a4_r, a4_i, a5_r, a5_i);
  B2_8(a6_r, a6_i, a7_r, a7_i);

  // Stage 2: Length 4
  B2_8(a0_r, a0_i, a2_r, a2_i);
  B2_8(a4_r, a4_i, a6_r, a6_i);

  #define ROT_I_8(xr, xi) { \
    __m256 tr = xr, ti = xi; \
    if constexpr (FORWARD) { \
      xr = ti; \
      xi = _mm256_sub_ps(v_zero, tr); \
    } else { \
      xr = _mm256_sub_ps(v_zero, ti); \
      xi = tr; \
    } \
  }

  ROT_I_8(a3_r, a3_i);
  ROT_I_8(a7_r, a7_i);

  B2_8(a1_r, a1_i, a3_r, a3_i);
  B2_8(a5_r, a5_i, a7_r, a7_i);

  // Stage 3: Length 8
  B2_8(a0_r, a0_i, a4_r, a4_i);

  // a5 * W8^1
  {
    __m256 ur = a5_r, ui = a5_i;
    if constexpr (FORWARD) {
      a5_r = _mm256_mul_ps(v_sqrt2_2, _mm256_add_ps(ur, ui));
      a5_i = _mm256_mul_ps(v_sqrt2_2, _mm256_sub_ps(ui, ur));
    } else {
      a5_r = _mm256_mul_ps(v_sqrt2_2, _mm256_sub_ps(ur, ui));
      a5_i = _mm256_mul_ps(v_sqrt2_2, _mm256_add_ps(ui, ur));
    }
  }
  B2_8(a1_r, a1_i, a5_r, a5_i);

  // a6 * W8^2
  ROT_I_8(a6_r, a6_i);
  B2_8(a2_r, a2_i, a6_r, a6_i);

  // a7 * W8^3
  {
    __m256 ur = a7_r, ui = a7_i;
    if constexpr (FORWARD) {
      a7_r = _mm256_mul_ps(v_sqrt2_2, _mm256_sub_ps(ui, ur));
      a7_i = _mm256_mul_ps(v_sqrt2_2, _mm256_sub_ps(_mm256_sub_ps(v_zero, ur), ui));
    } else {
      a7_r = _mm256_mul_ps(v_sqrt2_2, _mm256_sub_ps(_mm256_sub_ps(v_zero, ur), ui));
      a7_i = _mm256_mul_ps(v_sqrt2_2, _mm256_sub_ps(ur, ui));
    }
  }
  B2_8(a3_r, a3_i, a7_r, a7_i);

  #undef B2_8
  #undef ROT_I_8

  r[0] = a0_r; i[0] = a0_i;
  r[1] = a1_r; i[1] = a1_i;
  r[2] = a2_r; i[2] = a2_i;
  r[3] = a3_r; i[3] = a3_i;
  r[4] = a4_r; i[4] = a4_i;
  r[5] = a5_r; i[5] = a5_i;
  r[6] = a6_r; i[6] = a6_i;
  r[7] = a7_r; i[7] = a7_i;
}

} // namespace

void batch_fft8x8_r2c(std::size_t batch, const float* in, std::size_t in_dist,
                      std::size_t in_stride, cfloat* out, std::size_t out_dist,
                      std::size_t out_stride) noexcept {
  for (std::size_t b = 0; b < batch; ++b) {
    const float* blk_in = in + b * in_dist;
    cfloat* blk_out = out + b * out_dist;

    __m256 r[8];
    const __m256 vzero = _mm256_setzero_ps();
    for (int y = 0; y < 8; ++y) {
      r[y] = _mm256_loadu_ps(blk_in + y * in_stride);
    }

    // 1. Transpose: lane i becomes row i, register y becomes col y
    transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);

    // 2. Horizontal FFT on 8 rows
    __m256 i[8];
    for (int k = 0; k < 8; ++k) i[k] = vzero;
    fft8_8way_avx2<true>(r, i);

    // 3. Transpose back: register kx has 8 rows in lanes
    transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    transpose8x8_ps(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);

    // 4. Vertical FFT on 8 cols
    fft8_8way_avx2<true>(r, i);

    // 5. Store out 8x5 complex
    for (int ky = 0; ky < 8; ++ky) {
      alignas(32) float tr[8], ti[8];
      _mm256_store_ps(tr, r[ky]);
      _mm256_store_ps(ti, i[ky]);
      cfloat* dst = blk_out + ky * out_stride;
      for (int kx = 0; kx < 5; ++kx) {
        dst[kx] = cfloat(tr[kx], ti[kx]);
      }
    }
  }
}

void batch_fft8x8_c2r(std::size_t batch, const cfloat* in, std::size_t in_dist,
                      std::size_t in_stride, float* out, std::size_t out_dist,
                      std::size_t out_stride, float fct) noexcept {
  const __m256 vzero = _mm256_setzero_ps();
  const __m256 vfct = _mm256_set1_ps(fct);

  for (std::size_t b = 0; b < batch; ++b) {
    const cfloat* blk_in = in + b * in_dist;
    float* blk_out = out + b * out_dist;

    __m256 r[8], i[8];

    // 1. Load 8 rows of 5 complex columns directly into registers r[0..7]
    for (int ky = 0; ky < 8; ++ky) {
      const cfloat* row = blk_in + ky * in_stride;
      alignas(32) float tr[8] = {0}, ti[8] = {0};
      for (int kx = 0; kx < 5; ++kx) {
        tr[kx] = row[kx].real();
        ti[kx] = row[kx].imag();
      }
      r[ky] = _mm256_load_ps(tr);
      i[ky] = _mm256_load_ps(ti);
    }

    // 2. Vertical 1D IFFT along ky (the registers 0..7)
    fft8_8way_avx2<false>(r, i);

    // 3. Transpose so register index is kx, and lane is y:
    transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    transpose8x8_ps(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);

    // 4. Reconstruct Hermitian registers 5, 6, 7:
    r[5] = r[3]; i[5] = _mm256_sub_ps(vzero, i[3]);
    r[6] = r[2]; i[6] = _mm256_sub_ps(vzero, i[2]);
    r[7] = r[1]; i[7] = _mm256_sub_ps(vzero, i[1]);

    // 5. Horizontal 1D IFFT along kx (the registers 0..7)
    fft8_8way_avx2<false>(r, i);

    // 6. Scale by fct:
    for (int x = 0; x < 8; ++x) {
      r[x] = _mm256_mul_ps(r[x], vfct);
    }

    // 7. Transpose back so register index is y, and lane is x:
    transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);

    // 8. Store 8 rows:
    for (int y = 0; y < 8; ++y) {
      _mm256_storeu_ps(blk_out + y * out_stride, r[y]);
    }
  }
}

} // namespace neo_fft::codelet
#endif
