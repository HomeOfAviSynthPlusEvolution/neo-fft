#include "spectral/codelets/fft32x32.hpp"

#if defined(__AVX2__)
#include <immintrin.h>
#include <cmath>

namespace neo_fft::codelet {
namespace {

constexpr float PI = 3.14159265358979323846f;

static const float W16_R[8] = {
    std::cos(0.0f * PI / 8.0f), std::cos(1.0f * PI / 8.0f),
    std::cos(2.0f * PI / 8.0f), std::cos(3.0f * PI / 8.0f),
    std::cos(4.0f * PI / 8.0f), std::cos(5.0f * PI / 8.0f),
    std::cos(6.0f * PI / 8.0f), std::cos(7.0f * PI / 8.0f),
};

static const float W16_I[8] = {
    std::sin(0.0f * PI / 8.0f), std::sin(1.0f * PI / 8.0f),
    std::sin(2.0f * PI / 8.0f), std::sin(3.0f * PI / 8.0f),
    std::sin(4.0f * PI / 8.0f), std::sin(5.0f * PI / 8.0f),
    std::sin(6.0f * PI / 8.0f), std::sin(7.0f * PI / 8.0f),
};

static const float W32_R[16] = {
    std::cos(0.0f * PI / 16.0f),  std::cos(1.0f * PI / 16.0f),
    std::cos(2.0f * PI / 16.0f),  std::cos(3.0f * PI / 16.0f),
    std::cos(4.0f * PI / 16.0f),  std::cos(5.0f * PI / 16.0f),
    std::cos(6.0f * PI / 16.0f),  std::cos(7.0f * PI / 16.0f),
    std::cos(8.0f * PI / 16.0f),  std::cos(9.0f * PI / 16.0f),
    std::cos(10.0f * PI / 16.0f), std::cos(11.0f * PI / 16.0f),
    std::cos(12.0f * PI / 16.0f), std::cos(13.0f * PI / 16.0f),
    std::cos(14.0f * PI / 16.0f), std::cos(15.0f * PI / 16.0f),
};

static const float W32_I[16] = {
    std::sin(0.0f * PI / 16.0f),  std::sin(1.0f * PI / 16.0f),
    std::sin(2.0f * PI / 16.0f),  std::sin(3.0f * PI / 16.0f),
    std::sin(4.0f * PI / 16.0f),  std::sin(5.0f * PI / 16.0f),
    std::sin(6.0f * PI / 16.0f),  std::sin(7.0f * PI / 16.0f),
    std::sin(8.0f * PI / 16.0f),  std::sin(9.0f * PI / 16.0f),
    std::sin(10.0f * PI / 16.0f), std::sin(11.0f * PI / 16.0f),
    std::sin(12.0f * PI / 16.0f), std::sin(13.0f * PI / 16.0f),
    std::sin(14.0f * PI / 16.0f), std::sin(15.0f * PI / 16.0f),
};

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

#define CMUL_SCALAR(out_r, out_i, in_r, in_i, wr_val, wi_val, FORWARD) { \
  __m256 wr = _mm256_set1_ps(wr_val); \
  __m256 wi = _mm256_set1_ps((FORWARD) ? -(wi_val) : (wi_val)); \
  out_r = _mm256_fmsub_ps(in_r, wr, _mm256_mul_ps(in_i, wi)); \
  out_i = _mm256_fmadd_ps(in_r, wi, _mm256_mul_ps(in_i, wr)); \
}

#define MUL_J(out_r, out_i, in_r, in_i, FORWARD) { \
  if constexpr (FORWARD) { \
    out_r = in_i; \
    out_i = _mm256_sub_ps(_mm256_setzero_ps(), in_r); \
  } else { \
    out_r = _mm256_sub_ps(_mm256_setzero_ps(), in_i); \
    out_i = in_r; \
  } \
}

#define BUTTERFLY2(x_r, x_i, y_r, y_i) { \
  __m256 ur = x_r, ui = x_i; \
  __m256 vr = y_r, vi = y_i; \
  x_r = _mm256_add_ps(ur, vr); \
  x_i = _mm256_add_ps(ui, vi); \
  y_r = _mm256_sub_ps(ur, vr); \
  y_i = _mm256_sub_ps(ui, vi); \
}

template <bool FORWARD>
inline void fft16_8way_avx2(__m256* r, __m256* i) noexcept {
  __m256 a0_r = r[0],  a0_i = i[0];
  __m256 a1_r = r[8],  a1_i = i[8];
  __m256 a2_r = r[4],  a2_i = i[4];
  __m256 a3_r = r[12], a3_i = i[12];
  __m256 a4_r = r[2],  a4_i = i[2];
  __m256 a5_r = r[10], a5_i = i[10];
  __m256 a6_r = r[6],  a6_i = i[6];
  __m256 a7_r = r[14], a7_i = i[14];
  __m256 a8_r = r[1],  a8_i = i[1];
  __m256 a9_r = r[9],  a9_i = i[9];
  __m256 a10_r = r[5], a10_i = i[5];
  __m256 a11_r = r[13], a11_i = i[13];
  __m256 a12_r = r[3], a12_i = i[3];
  __m256 a13_r = r[11], a13_i = i[11];
  __m256 a14_r = r[7], a14_i = i[7];
  __m256 a15_r = r[15], a15_i = i[15];

  BUTTERFLY2(a0_r, a0_i, a1_r, a1_i);
  BUTTERFLY2(a2_r, a2_i, a3_r, a3_i);
  BUTTERFLY2(a4_r, a4_i, a5_r, a5_i);
  BUTTERFLY2(a6_r, a6_i, a7_r, a7_i);
  BUTTERFLY2(a8_r, a8_i, a9_r, a9_i);
  BUTTERFLY2(a10_r, a10_i, a11_r, a11_i);
  BUTTERFLY2(a12_r, a12_i, a13_r, a13_i);
  BUTTERFLY2(a14_r, a14_i, a15_r, a15_i);

  #define BUTTERFLY4(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i) { \
    BUTTERFLY2(x0_r, x0_i, x2_r, x2_i); \
    __m256 v1r, v1i; \
    MUL_J(v1r, v1i, x3_r, x3_i, FORWARD); \
    __m256 u1r = x1_r, u1i = x1_i; \
    x1_r = _mm256_add_ps(u1r, v1r); \
    x1_i = _mm256_add_ps(u1i, v1i); \
    x3_r = _mm256_sub_ps(u1r, v1r); \
    x3_i = _mm256_sub_ps(u1i, v1i); \
  }

  BUTTERFLY4(a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i);
  BUTTERFLY4(a4_r, a4_i, a5_r, a5_i, a6_r, a6_i, a7_r, a7_i);
  BUTTERFLY4(a8_r, a8_i, a9_r, a9_i, a10_r, a10_i, a11_r, a11_i);
  BUTTERFLY4(a12_r, a12_i, a13_r, a13_i, a14_r, a14_i, a15_r, a15_i);
  #undef BUTTERFLY4

  #define BUTTERFLY8(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, \
                     x4_r, x4_i, x5_r, x5_i, x6_r, x6_i, x7_r, x7_i) { \
    BUTTERFLY2(x0_r, x0_i, x4_r, x4_i); \
    __m256 v1r, v1i; \
    CMUL_SCALAR(v1r, v1i, x5_r, x5_i, W16_R[2], W16_I[2], FORWARD); \
    __m256 u1r = x1_r, u1i = x1_i; \
    x1_r = _mm256_add_ps(u1r, v1r); \
    x1_i = _mm256_add_ps(u1i, v1i); \
    x5_r = _mm256_sub_ps(u1r, v1r); \
    x5_i = _mm256_sub_ps(u1i, v1i); \
    __m256 v2r, v2i; \
    MUL_J(v2r, v2i, x6_r, x6_i, FORWARD); \
    __m256 u2r = x2_r, u2i = x2_i; \
    x2_r = _mm256_add_ps(u2r, v2r); \
    x2_i = _mm256_add_ps(u2i, v2i); \
    x6_r = _mm256_sub_ps(u2r, v2r); \
    x6_i = _mm256_sub_ps(u2i, v2i); \
    __m256 v3r, v3i; \
    CMUL_SCALAR(v3r, v3i, x7_r, x7_i, W16_R[6], W16_I[6], FORWARD); \
    __m256 u3r = x3_r, u3i = x3_i; \
    x3_r = _mm256_add_ps(u3r, v3r); \
    x3_i = _mm256_add_ps(u3i, v3i); \
    x7_r = _mm256_sub_ps(u3r, v3r); \
    x7_i = _mm256_sub_ps(u3i, v3i); \
  }

  BUTTERFLY8(a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i,
             a4_r, a4_i, a5_r, a5_i, a6_r, a6_i, a7_r, a7_i);
  BUTTERFLY8(a8_r, a8_i, a9_r, a9_i, a10_r, a10_i, a11_r, a11_i,
             a12_r, a12_i, a13_r, a13_i, a14_r, a14_i, a15_r, a15_i);
  #undef BUTTERFLY8

  BUTTERFLY2(a0_r, a0_i, a8_r, a8_i);
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a9_r, a9_i, W16_R[1], W16_I[1], FORWARD);
    __m256 ur = a1_r, ui = a1_i;
    a1_r = _mm256_add_ps(ur, vr); a1_i = _mm256_add_ps(ui, vi);
    a9_r = _mm256_sub_ps(ur, vr); a9_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a10_r, a10_i, W16_R[2], W16_I[2], FORWARD);
    __m256 ur = a2_r, ui = a2_i;
    a2_r = _mm256_add_ps(ur, vr); a2_i = _mm256_add_ps(ui, vi);
    a10_r = _mm256_sub_ps(ur, vr); a10_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a11_r, a11_i, W16_R[3], W16_I[3], FORWARD);
    __m256 ur = a3_r, ui = a3_i;
    a3_r = _mm256_add_ps(ur, vr); a3_i = _mm256_add_ps(ui, vi);
    a11_r = _mm256_sub_ps(ur, vr); a11_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    MUL_J(vr, vi, a12_r, a12_i, FORWARD);
    __m256 ur = a4_r, ui = a4_i;
    a4_r = _mm256_add_ps(ur, vr); a4_i = _mm256_add_ps(ui, vi);
    a12_r = _mm256_sub_ps(ur, vr); a12_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a13_r, a13_i, W16_R[5], W16_I[5], FORWARD);
    __m256 ur = a5_r, ui = a5_i;
    a5_r = _mm256_add_ps(ur, vr); a5_i = _mm256_add_ps(ui, vi);
    a13_r = _mm256_sub_ps(ur, vr); a13_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a14_r, a14_i, W16_R[6], W16_I[6], FORWARD);
    __m256 ur = a6_r, ui = a6_i;
    a6_r = _mm256_add_ps(ur, vr); a6_i = _mm256_add_ps(ui, vi);
    a14_r = _mm256_sub_ps(ur, vr); a14_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a15_r, a15_i, W16_R[7], W16_I[7], FORWARD);
    __m256 ur = a7_r, ui = a7_i;
    a7_r = _mm256_add_ps(ur, vr); a7_i = _mm256_add_ps(ui, vi);
    a15_r = _mm256_sub_ps(ur, vr); a15_i = _mm256_sub_ps(ui, vi);
  }

  r[0] = a0_r; i[0] = a0_i;
  r[1] = a1_r; i[1] = a1_i;
  r[2] = a2_r; i[2] = a2_i;
  r[3] = a3_r; i[3] = a3_i;
  r[4] = a4_r; i[4] = a4_i;
  r[5] = a5_r; i[5] = a5_i;
  r[6] = a6_r; i[6] = a6_i;
  r[7] = a7_r; i[7] = a7_i;
  r[8] = a8_r; i[8] = a8_i;
  r[9] = a9_r; i[9] = a9_i;
  r[10] = a10_r; i[10] = a10_i;
  r[11] = a11_r; i[11] = a11_i;
  r[12] = a12_r; i[12] = a12_i;
  r[13] = a13_r; i[13] = a13_i;
  r[14] = a14_r; i[14] = a14_i;
  r[15] = a15_r; i[15] = a15_i;
}

template <bool FORWARD>
inline void fft32_8way_avx2(__m256* r, __m256* i) noexcept {
  __m256 er[16], ei[16], or_[16], oi[16];
  for (int m = 0; m < 16; ++m) {
    er[m] = r[2 * m];
    ei[m] = i[2 * m];
    or_[m] = r[2 * m + 1];
    oi[m] = i[2 * m + 1];
  }

  fft16_8way_avx2<FORWARD>(er, ei);
  fft16_8way_avx2<FORWARD>(or_, oi);

  for (int k = 0; k < 16; ++k) {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, or_[k], oi[k], W32_R[k], W32_I[k], FORWARD);
    r[k]      = _mm256_add_ps(er[k], vr);
    i[k]      = _mm256_add_ps(ei[k], vi);
    r[k + 16] = _mm256_sub_ps(er[k], vr);
    i[k + 16] = _mm256_sub_ps(ei[k], vi);
  }
}

#undef BUTTERFLY2
#undef CMUL_SCALAR
#undef MUL_J

inline void batch_r2c_avx2(std::size_t batch, const float* in, std::size_t in_dist,
                           std::size_t in_stride, cfloat* out, std::size_t out_dist,
                           std::size_t out_stride) noexcept {
  for (std::size_t b = 0; b < batch; ++b) {
    const float* blk_in = in + b * in_dist;
    cfloat* blk_out = out + b * out_dist;

    alignas(32) float inter_r[32][32];
    alignas(32) float inter_i[32][32];

    // 1. Horizontal 1D FFT on 32 rows (4 groups of 8 rows)
    for (int g = 0; g < 4; ++g) {
      const int y_base = g * 8;
      __m256 r[32], i[32];
      const __m256 vzero = _mm256_setzero_ps();

      for (int chunk = 0; chunk < 4; ++chunk) {
        const int x_base = chunk * 8;
        __m256 b0 = _mm256_loadu_ps(blk_in + (y_base + 0) * in_stride + x_base);
        __m256 b1 = _mm256_loadu_ps(blk_in + (y_base + 1) * in_stride + x_base);
        __m256 b2 = _mm256_loadu_ps(blk_in + (y_base + 2) * in_stride + x_base);
        __m256 b3 = _mm256_loadu_ps(blk_in + (y_base + 3) * in_stride + x_base);
        __m256 b4 = _mm256_loadu_ps(blk_in + (y_base + 4) * in_stride + x_base);
        __m256 b5 = _mm256_loadu_ps(blk_in + (y_base + 5) * in_stride + x_base);
        __m256 b6 = _mm256_loadu_ps(blk_in + (y_base + 6) * in_stride + x_base);
        __m256 b7 = _mm256_loadu_ps(blk_in + (y_base + 7) * in_stride + x_base);

        transpose8x8_ps(b0, b1, b2, b3, b4, b5, b6, b7);

        r[x_base + 0] = b0; r[x_base + 1] = b1; r[x_base + 2] = b2; r[x_base + 3] = b3;
        r[x_base + 4] = b4; r[x_base + 5] = b5; r[x_base + 6] = b6; r[x_base + 7] = b7;
      }
      for (int k = 0; k < 32; ++k) {
        i[k] = vzero;
      }

      fft32_8way_avx2<true>(r, i);

      // In-register transpose bins 0..7 and 8..15
      transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
      transpose8x8_ps(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);
      transpose8x8_ps(r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15]);
      transpose8x8_ps(i[8], i[9], i[10], i[11], i[12], i[13], i[14], i[15]);

      alignas(32) float tr16[8];
      _mm256_store_ps(tr16, r[16]);

      for (int y = 0; y < 8; ++y) {
        _mm256_store_ps(&inter_r[y_base + y][0], r[y]);
        _mm256_store_ps(&inter_i[y_base + y][0], i[y]);
        _mm256_store_ps(&inter_r[y_base + y][8], r[y + 8]);
        _mm256_store_ps(&inter_i[y_base + y][8], i[y + 8]);
        inter_r[y_base + y][16] = tr16[y];
        inter_i[y_base + y][16] = 0.0f;
      }
    }

    // 2. Vertical 1D FFT on 17 columns
    // Columns 0..7 (8-way parallel)
    {
      __m256 r[32], i[32];
      for (int y = 0; y < 32; ++y) {
        r[y] = _mm256_load_ps(&inter_r[y][0]);
        i[y] = _mm256_load_ps(&inter_i[y][0]);
      }
      fft32_8way_avx2<true>(r, i);
      for (int ky = 0; ky < 32; ++ky) {
        cfloat* dst = blk_out + ky * out_stride;
        __m256 xy_lo = _mm256_unpacklo_ps(r[ky], i[ky]);
        __m256 xy_hi = _mm256_unpackhi_ps(r[ky], i[ky]);
        __m256 c0_3  = _mm256_permute2f128_ps(xy_lo, xy_hi, 0x20);
        __m256 c4_7  = _mm256_permute2f128_ps(xy_lo, xy_hi, 0x31);
        _mm256_storeu_ps(reinterpret_cast<float*>(dst), c0_3);
        _mm256_storeu_ps(reinterpret_cast<float*>(dst + 4), c4_7);
      }
    }
    // Columns 8..15 (8-way parallel)
    {
      __m256 r[32], i[32];
      for (int y = 0; y < 32; ++y) {
        r[y] = _mm256_load_ps(&inter_r[y][8]);
        i[y] = _mm256_load_ps(&inter_i[y][8]);
      }
      fft32_8way_avx2<true>(r, i);
      for (int ky = 0; ky < 32; ++ky) {
        cfloat* dst = blk_out + ky * out_stride + 8;
        __m256 xy_lo = _mm256_unpacklo_ps(r[ky], i[ky]);
        __m256 xy_hi = _mm256_unpackhi_ps(r[ky], i[ky]);
        __m256 c0_3  = _mm256_permute2f128_ps(xy_lo, xy_hi, 0x20);
        __m256 c4_7  = _mm256_permute2f128_ps(xy_lo, xy_hi, 0x31);
        _mm256_storeu_ps(reinterpret_cast<float*>(dst), c0_3);
        _mm256_storeu_ps(reinterpret_cast<float*>(dst + 4), c4_7);
      }
    }
    // Column 16 (single column)
    {
      __m256 r[32], i[32];
      for (int y = 0; y < 32; ++y) {
        r[y] = _mm256_set1_ps(inter_r[y][16]);
        i[y] = _mm256_set1_ps(inter_i[y][16]);
      }
      fft32_8way_avx2<true>(r, i);
      for (int ky = 0; ky < 32; ++ky) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[ky]);
        _mm256_store_ps(ti, i[ky]);
        blk_out[ky * out_stride + 16] = cfloat(tr[0], ti[0]);
      }
    }
  }
}

inline void batch_c2r_avx2(std::size_t batch, const cfloat* in, std::size_t in_dist,
                           std::size_t in_stride, float* out, std::size_t out_dist,
                           std::size_t out_stride, float fct) noexcept {
  const __m256 vfct = _mm256_set1_ps(fct);

  for (std::size_t b = 0; b < batch; ++b) {
    const cfloat* blk_in = in + b * in_dist;
    float* blk_out = out + b * out_dist;

    alignas(32) float inter_r[32][32];
    alignas(32) float inter_i[32][32];

    // 1. Vertical 1D IFFT on 17 columns
    // Columns 0..7
    {
      __m256 r[32], i[32];
      for (int y = 0; y < 32; ++y) {
        const cfloat* src = blk_in + y * in_stride;
        __m256 c0_3 = _mm256_loadu_ps(reinterpret_cast<const float*>(src));
        __m256 c4_7 = _mm256_loadu_ps(reinterpret_cast<const float*>(src + 4));
        __m256 s_lo = _mm256_shuffle_ps(c0_3, c4_7, _MM_SHUFFLE(2, 0, 2, 0));
        __m256 s_hi = _mm256_shuffle_ps(c0_3, c4_7, _MM_SHUFFLE(3, 1, 3, 1));
        r[y] = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(s_lo), _MM_SHUFFLE(3, 1, 2, 0)));
        i[y] = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(s_hi), _MM_SHUFFLE(3, 1, 2, 0)));
      }
      fft32_8way_avx2<false>(r, i);
      for (int ky = 0; ky < 32; ++ky) {
        _mm256_store_ps(&inter_r[ky][0], r[ky]);
        _mm256_store_ps(&inter_i[ky][0], i[ky]);
      }
    }
    // Columns 8..15
    {
      __m256 r[32], i[32];
      for (int y = 0; y < 32; ++y) {
        const cfloat* src = blk_in + y * in_stride + 8;
        __m256 c0_3 = _mm256_loadu_ps(reinterpret_cast<const float*>(src));
        __m256 c4_7 = _mm256_loadu_ps(reinterpret_cast<const float*>(src + 4));
        __m256 s_lo = _mm256_shuffle_ps(c0_3, c4_7, _MM_SHUFFLE(2, 0, 2, 0));
        __m256 s_hi = _mm256_shuffle_ps(c0_3, c4_7, _MM_SHUFFLE(3, 1, 3, 1));
        r[y] = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(s_lo), _MM_SHUFFLE(3, 1, 2, 0)));
        i[y] = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(s_hi), _MM_SHUFFLE(3, 1, 2, 0)));
      }
      fft32_8way_avx2<false>(r, i);
      for (int ky = 0; ky < 32; ++ky) {
        _mm256_store_ps(&inter_r[ky][8], r[ky]);
        _mm256_store_ps(&inter_i[ky][8], i[ky]);
      }
    }
    // Column 16
    {
      __m256 r[32], i[32];
      for (int y = 0; y < 32; ++y) {
        const cfloat* src = blk_in + y * in_stride;
        r[y] = _mm256_set1_ps(src[16].real());
        i[y] = _mm256_set1_ps(src[16].imag());
      }
      fft32_8way_avx2<false>(r, i);
      for (int ky = 0; ky < 32; ++ky) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[ky]);
        _mm256_store_ps(ti, i[ky]);
        inter_r[ky][16] = tr[0];
        inter_i[ky][16] = ti[0];
      }
    }

    // 2. Horizontal 1D IFFT on 32 rows directly to blk_out
    for (int g = 0; g < 4; ++g) {
      const int y_base = g * 8;
      __m256 r[32], i[32];

      // Load bins 0..7
      for (int y = 0; y < 8; ++y) {
        r[y] = _mm256_load_ps(&inter_r[y_base + y][0]);
        i[y] = _mm256_load_ps(&inter_i[y_base + y][0]);
      }
      transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
      transpose8x8_ps(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);

      // Load bins 8..15
      for (int y = 0; y < 8; ++y) {
        r[8 + y] = _mm256_load_ps(&inter_r[y_base + y][8]);
        i[8 + y] = _mm256_load_ps(&inter_i[y_base + y][8]);
      }
      transpose8x8_ps(r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15]);
      transpose8x8_ps(i[8], i[9], i[10], i[11], i[12], i[13], i[14], i[15]);

      // Load bin 16
      alignas(32) float tr16[8], ti16[8];
      for (int y = 0; y < 8; ++y) {
        tr16[y] = inter_r[y_base + y][16];
        ti16[y] = inter_i[y_base + y][16];
      }
      r[16] = _mm256_load_ps(tr16);
      i[16] = _mm256_load_ps(ti16);

      // Reconstruct bins 17..31 with conjugate symmetry
      const __m256 vzero = _mm256_setzero_ps();
      for (int kx = 1; kx < 16; ++kx) {
        r[32 - kx] = r[kx];
        i[32 - kx] = _mm256_sub_ps(vzero, i[kx]);
      }

      fft32_8way_avx2<false>(r, i);

      for (int k = 0; k < 32; ++k) {
        r[k] = _mm256_mul_ps(r[k], vfct);
      }

      transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
      transpose8x8_ps(r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15]);
      transpose8x8_ps(r[16], r[17], r[18], r[19], r[20], r[21], r[22], r[23]);
      transpose8x8_ps(r[24], r[25], r[26], r[27], r[28], r[29], r[30], r[31]);

      for (int y = 0; y < 8; ++y) {
        float* dst = blk_out + (y_base + y) * out_stride;
        _mm256_storeu_ps(dst + 0,  r[y]);
        _mm256_storeu_ps(dst + 8,  r[y + 8]);
        _mm256_storeu_ps(dst + 16, r[y + 16]);
        _mm256_storeu_ps(dst + 24, r[y + 24]);
      }
    }
  }
}

} // namespace

void batch_fft32x32_r2c(std::size_t batch, const float* in, std::size_t in_dist,
                        std::size_t in_stride, cfloat* out, std::size_t out_dist,
                        std::size_t out_stride) noexcept {
  batch_r2c_avx2(batch, in, in_dist, in_stride, out, out_dist, out_stride);
}

void batch_fft32x32_c2r(std::size_t batch, const cfloat* in, std::size_t in_dist,
                        std::size_t in_stride, float* out, std::size_t out_dist,
                        std::size_t out_stride, float fct) noexcept {
  batch_c2r_avx2(batch, in, in_dist, in_stride, out, out_dist, out_stride, fct);
}

} // namespace neo_fft::codelet
#endif
