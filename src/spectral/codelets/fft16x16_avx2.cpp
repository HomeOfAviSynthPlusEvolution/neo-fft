#include "spectral/codelets/fft16x16.hpp"

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

  #define BUTTERFLY2(x_r, x_i, y_r, y_i) { \
    __m256 ur = x_r, ui = x_i; \
    __m256 vr = y_r, vi = y_i; \
    x_r = _mm256_add_ps(ur, vr); \
    x_i = _mm256_add_ps(ui, vi); \
    y_r = _mm256_sub_ps(ur, vr); \
    y_i = _mm256_sub_ps(ui, vi); \
  }

  BUTTERFLY2(a0_r, a0_i, a1_r, a1_i);
  BUTTERFLY2(a2_r, a2_i, a3_r, a3_i);
  BUTTERFLY2(a4_r, a4_i, a5_r, a5_i);
  BUTTERFLY2(a6_r, a6_i, a7_r, a7_i);
  BUTTERFLY2(a8_r, a8_i, a9_r, a9_i);
  BUTTERFLY2(a10_r, a10_i, a11_r, a11_i);
  BUTTERFLY2(a12_r, a12_i, a13_r, a13_i);
  BUTTERFLY2(a14_r, a14_i, a15_r, a15_i);

  #define CMUL_SCALAR(out_r, out_i, in_r, in_i, wr_val, wi_val) { \
    __m256 wr = _mm256_set1_ps(wr_val); \
    __m256 wi = _mm256_set1_ps(FORWARD ? -(wi_val) : (wi_val)); \
    out_r = _mm256_fmsub_ps(in_r, wr, _mm256_mul_ps(in_i, wi)); \
    out_i = _mm256_fmadd_ps(in_r, wi, _mm256_mul_ps(in_i, wr)); \
  }

  #define MUL_J(out_r, out_i, in_r, in_i) { \
    if constexpr (FORWARD) { \
      out_r = in_i; \
      out_i = _mm256_sub_ps(_mm256_setzero_ps(), in_r); \
    } else { \
      out_r = _mm256_sub_ps(_mm256_setzero_ps(), in_i); \
      out_i = in_r; \
    } \
  }

  #define BUTTERFLY4(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i) { \
    BUTTERFLY2(x0_r, x0_i, x2_r, x2_i); \
    __m256 v1r, v1i; \
    MUL_J(v1r, v1i, x3_r, x3_i); \
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

  #define BUTTERFLY8(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, \
                     x4_r, x4_i, x5_r, x5_i, x6_r, x6_i, x7_r, x7_i) { \
    BUTTERFLY2(x0_r, x0_i, x4_r, x4_i); \
    __m256 v1r, v1i; \
    CMUL_SCALAR(v1r, v1i, x5_r, x5_i, W16_R[2], W16_I[2]); \
    __m256 u1r = x1_r, u1i = x1_i; \
    x1_r = _mm256_add_ps(u1r, v1r); \
    x1_i = _mm256_add_ps(u1i, v1i); \
    x5_r = _mm256_sub_ps(u1r, v1r); \
    x5_i = _mm256_sub_ps(u1i, v1i); \
    __m256 v2r, v2i; \
    MUL_J(v2r, v2i, x6_r, x6_i); \
    __m256 u2r = x2_r, u2i = x2_i; \
    x2_r = _mm256_add_ps(u2r, v2r); \
    x2_i = _mm256_add_ps(u2i, v2i); \
    x6_r = _mm256_sub_ps(u2r, v2r); \
    x6_i = _mm256_sub_ps(u2i, v2i); \
    __m256 v3r, v3i; \
    CMUL_SCALAR(v3r, v3i, x7_r, x7_i, W16_R[6], W16_I[6]); \
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

  BUTTERFLY2(a0_r, a0_i, a8_r, a8_i);
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a9_r, a9_i, W16_R[1], W16_I[1]);
    __m256 ur = a1_r, ui = a1_i;
    a1_r = _mm256_add_ps(ur, vr); a1_i = _mm256_add_ps(ui, vi);
    a9_r = _mm256_sub_ps(ur, vr); a9_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a10_r, a10_i, W16_R[2], W16_I[2]);
    __m256 ur = a2_r, ui = a2_i;
    a2_r = _mm256_add_ps(ur, vr); a2_i = _mm256_add_ps(ui, vi);
    a10_r = _mm256_sub_ps(ur, vr); a10_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a11_r, a11_i, W16_R[3], W16_I[3]);
    __m256 ur = a3_r, ui = a3_i;
    a3_r = _mm256_add_ps(ur, vr); a3_i = _mm256_add_ps(ui, vi);
    a11_r = _mm256_sub_ps(ur, vr); a11_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    MUL_J(vr, vi, a12_r, a12_i);
    __m256 ur = a4_r, ui = a4_i;
    a4_r = _mm256_add_ps(ur, vr); a4_i = _mm256_add_ps(ui, vi);
    a12_r = _mm256_sub_ps(ur, vr); a12_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a13_r, a13_i, W16_R[5], W16_I[5]);
    __m256 ur = a5_r, ui = a5_i;
    a5_r = _mm256_add_ps(ur, vr); a5_i = _mm256_add_ps(ui, vi);
    a13_r = _mm256_sub_ps(ur, vr); a13_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a14_r, a14_i, W16_R[6], W16_I[6]);
    __m256 ur = a6_r, ui = a6_i;
    a6_r = _mm256_add_ps(ur, vr); a6_i = _mm256_add_ps(ui, vi);
    a14_r = _mm256_sub_ps(ur, vr); a14_i = _mm256_sub_ps(ui, vi);
  }
  {
    __m256 vr, vi;
    CMUL_SCALAR(vr, vi, a15_r, a15_i, W16_R[7], W16_I[7]);
    __m256 ur = a7_r, ui = a7_i;
    a7_r = _mm256_add_ps(ur, vr); a7_i = _mm256_add_ps(ui, vi);
    a15_r = _mm256_sub_ps(ur, vr); a15_i = _mm256_sub_ps(ui, vi);
  }

  #undef BUTTERFLY2
  #undef BUTTERFLY4
  #undef BUTTERFLY8
  #undef CMUL_SCALAR
  #undef MUL_J

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

inline void batch_r2c_avx2(std::size_t batch, const float* in, std::size_t in_dist,
                           std::size_t in_stride, cfloat* out, std::size_t out_dist,
                           std::size_t out_stride) noexcept {
  for (std::size_t b = 0; b < batch; ++b) {
    const float* blk_in = in + b * in_dist;
    cfloat* blk_out = out + b * out_dist;

    alignas(32) float inter_r[16][16];
    alignas(32) float inter_i[16][16];

    // 1. Horizontal 1D FFT on 16 rows (2 groups of 8 rows)
    for (int g = 0; g < 2; ++g) {
      const int y_base = g * 8;
      __m256 r0 = _mm256_loadu_ps(blk_in + (y_base + 0) * in_stride);
      __m256 r1 = _mm256_loadu_ps(blk_in + (y_base + 1) * in_stride);
      __m256 r2 = _mm256_loadu_ps(blk_in + (y_base + 2) * in_stride);
      __m256 r3 = _mm256_loadu_ps(blk_in + (y_base + 3) * in_stride);
      __m256 r4 = _mm256_loadu_ps(blk_in + (y_base + 4) * in_stride);
      __m256 r5 = _mm256_loadu_ps(blk_in + (y_base + 5) * in_stride);
      __m256 r6 = _mm256_loadu_ps(blk_in + (y_base + 6) * in_stride);
      __m256 r7 = _mm256_loadu_ps(blk_in + (y_base + 7) * in_stride);

      __m256 r8  = _mm256_loadu_ps(blk_in + (y_base + 0) * in_stride + 8);
      __m256 r9  = _mm256_loadu_ps(blk_in + (y_base + 1) * in_stride + 8);
      __m256 r10 = _mm256_loadu_ps(blk_in + (y_base + 2) * in_stride + 8);
      __m256 r11 = _mm256_loadu_ps(blk_in + (y_base + 3) * in_stride + 8);
      __m256 r12 = _mm256_loadu_ps(blk_in + (y_base + 4) * in_stride + 8);
      __m256 r13 = _mm256_loadu_ps(blk_in + (y_base + 5) * in_stride + 8);
      __m256 r14 = _mm256_loadu_ps(blk_in + (y_base + 6) * in_stride + 8);
      __m256 r15 = _mm256_loadu_ps(blk_in + (y_base + 7) * in_stride + 8);

      transpose8x8_ps(r0, r1, r2, r3, r4, r5, r6, r7);
      transpose8x8_ps(r8, r9, r10, r11, r12, r13, r14, r15);

      __m256 r[16] = {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15};
      __m256 i[16];
      const __m256 vzero = _mm256_setzero_ps();
      for (int k = 0; k < 16; ++k) {
        i[k] = vzero;
      }

      fft16_8way_avx2<true>(r, i);

      for (int kx = 0; kx < 9; ++kx) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[kx]);
        _mm256_store_ps(ti, i[kx]);
        for (int y = 0; y < 8; ++y) {
          inter_r[y_base + y][kx] = tr[y];
          inter_i[y_base + y][kx] = ti[y];
        }
      }
    }

    // 2. Vertical 1D FFT on 9 columns
    // Columns 0..7 (8-way parallel)
    {
      __m256 r[16], i[16];
      for (int y = 0; y < 16; ++y) {
        r[y] = _mm256_load_ps(inter_r[y]);
        i[y] = _mm256_load_ps(inter_i[y]);
      }
      fft16_8way_avx2<true>(r, i);
      for (int ky = 0; ky < 16; ++ky) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[ky]);
        _mm256_store_ps(ti, i[ky]);
        cfloat* dst = blk_out + ky * out_stride;
        for (int kx = 0; kx < 8; ++kx) {
          dst[kx] = cfloat(tr[kx], ti[kx]);
        }
      }
    }
    // Column 8 (single column)
    {
      __m256 r[16], i[16];
      for (int y = 0; y < 16; ++y) {
        r[y] = _mm256_set1_ps(inter_r[y][8]);
        i[y] = _mm256_set1_ps(inter_i[y][8]);
      }
      fft16_8way_avx2<true>(r, i);
      for (int ky = 0; ky < 16; ++ky) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[ky]);
        _mm256_store_ps(ti, i[ky]);
        cfloat* dst = blk_out + ky * out_stride;
        dst[8] = cfloat(tr[0], ti[0]);
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

    alignas(32) float inter_r[16][16];
    alignas(32) float inter_i[16][16];

    // 1. Vertical 1D IFFT on 9 columns
    {
      __m256 r[16], i[16];
      for (int y = 0; y < 16; ++y) {
        const cfloat* src = blk_in + y * in_stride;
        alignas(32) float tr[8], ti[8];
        for (int kx = 0; kx < 8; ++kx) {
          tr[kx] = src[kx].real();
          ti[kx] = src[kx].imag();
        }
        r[y] = _mm256_load_ps(tr);
        i[y] = _mm256_load_ps(ti);
      }
      fft16_8way_avx2<false>(r, i);
      for (int ky = 0; ky < 16; ++ky) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[ky]);
        _mm256_store_ps(ti, i[ky]);
        for (int kx = 0; kx < 8; ++kx) {
          inter_r[ky][kx] = tr[kx];
          inter_i[ky][kx] = ti[kx];
        }
      }
    }
    {
      __m256 r[16], i[16];
      for (int y = 0; y < 16; ++y) {
        const cfloat* src = blk_in + y * in_stride;
        r[y] = _mm256_set1_ps(src[8].real());
        i[y] = _mm256_set1_ps(src[8].imag());
      }
      fft16_8way_avx2<false>(r, i);
      for (int ky = 0; ky < 16; ++ky) {
        alignas(32) float tr[8], ti[8];
        _mm256_store_ps(tr, r[ky]);
        _mm256_store_ps(ti, i[ky]);
        inter_r[ky][8] = tr[0];
        inter_i[ky][8] = ti[0];
      }
    }

    // 2. Horizontal 1D IFFT directly to blk_out
    for (int g = 0; g < 2; ++g) {
      const int y_base = g * 8;
      __m256 r[16], i[16];
      for (int kx = 0; kx < 9; ++kx) {
        alignas(32) float tr[8], ti[8];
        for (int y = 0; y < 8; ++y) {
          tr[y] = inter_r[y_base + y][kx];
          ti[y] = inter_i[y_base + y][kx];
        }
        r[kx] = _mm256_load_ps(tr);
        i[kx] = _mm256_load_ps(ti);
      }
      const __m256 vzero = _mm256_setzero_ps();
      for (int kx = 1; kx < 8; ++kx) {
        r[16 - kx] = r[kx];
        i[16 - kx] = _mm256_sub_ps(vzero, i[kx]);
      }

      fft16_8way_avx2<false>(r, i);

      for (int x = 0; x < 16; ++x) {
        r[x] = _mm256_mul_ps(r[x], vfct);
      }

      transpose8x8_ps(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
      transpose8x8_ps(r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15]);

      for (int y = 0; y < 8; ++y) {
        float* dst = blk_out + (y_base + y) * out_stride;
        _mm256_storeu_ps(dst, r[y]);
        _mm256_storeu_ps(dst + 8, r[y + 8]);
      }
    }
  }
}

} // namespace

void batch_fft16x16_r2c(std::size_t batch, const float* in, std::size_t in_dist,
                        std::size_t in_stride, cfloat* out, std::size_t out_dist,
                        std::size_t out_stride) noexcept {
  batch_r2c_avx2(batch, in, in_dist, in_stride, out, out_dist, out_stride);
}

void batch_fft16x16_c2r(std::size_t batch, const cfloat* in, std::size_t in_dist,
                        std::size_t in_stride, float* out, std::size_t out_dist,
                        std::size_t out_stride, float fct) noexcept {
  batch_c2r_avx2(batch, in, in_dist, in_stride, out, out_dist, out_stride, fct);
}

} // namespace neo_fft::codelet
#endif
