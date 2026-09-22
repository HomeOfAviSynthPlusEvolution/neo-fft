#include "kernels/spectral.hpp"
#include "kernels/spatial.hpp"
#include "../test.hpp"
#include <algorithm>
#include <cstring>
#include <memory>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
using namespace neo_fft;
struct Guarded {
  std::complex<float>* data = nullptr;
  void* allocation = nullptr;
  std::vector<std::complex<float>> fallback;
  explicit Guarded(std::size_t count) {
#if defined(_WIN32)
    SYSTEM_INFO info{};
    GetSystemInfo(&info);
    CHECK(count * sizeof(std::complex<float>) <= info.dwPageSize);
    allocation = VirtualAlloc(nullptr, 2 * info.dwPageSize, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    CHECK(allocation);
    auto* boundary = static_cast<unsigned char*>(allocation) + info.dwPageSize;
    DWORD previous = 0;
    CHECK(VirtualProtect(boundary, info.dwPageSize, PAGE_NOACCESS, &previous));
    data = reinterpret_cast<std::complex<float>*>(boundary - count * sizeof(std::complex<float>));
    for (std::size_t i = 0; i < count; ++i)
      new (data + i) std::complex<float>();
#else
    fallback.resize(count);
    data = fallback.data();
#endif
  }
  ~Guarded() {
#if defined(_WIN32)
    if (allocation)
      VirtualFree(allocation, 0, MEM_RELEASE);
#endif
  }
  Guarded(const Guarded&) = delete;
  Guarded& operator=(const Guarded&) = delete;
};
int main() {
  try {
    // Every offset, vector boundary and odd tail ends exactly at a guard page.
    // Preserve input bits, reject +/-Inf and NaN even away from the first lane.
    for (int opt : {0, 1}) {
      const auto check = select_spatial(opt).validate_finite;
      check(nullptr, 0);
      for (std::size_t count = 1; count <= 257; ++count) {
        Guarded storage(count);
        auto* values = reinterpret_cast<float*>(storage.data) + count;
        for (std::size_t i = 0; i < count; ++i)
          values[i] = i % 2 ? -std::numeric_limits<float>::max() : std::numeric_limits<float>::denorm_min();
        const std::vector<float> original(values, values + count);
        check(values, count);
        CHECK(std::memcmp(values, original.data(), count * sizeof(float)) == 0);
        for (std::size_t i = 0; i < count; ++i) {
          for (float invalid : {INFINITY, -INFINITY, NAN}) {
            values[i] = invalid;
            bool failed = false;
            try { check(values, count); }
            catch (const std::runtime_error& e) {
              CHECK(std::string(e.what()) == "non-finite sample or intermediate");
              failed = true;
            }
            CHECK(failed);
          }
          values[i] = original[i];
        }
      }
    }
    for(int opt:{0,1}) {
      std::vector<std::complex<float>> values(32,{1,0});
      SpectralParams p;p.type=4;p.a=1;p.low=p.high=1e20f;
      rejects([&]{select_spectral(opt)(values.data(),nullptr,values.size(),0,p);});
    }
    // Standalone kernels retain shape admission; the Plan-only path skips
    // shapes, not numeric overflow checks.
    for(int opt:{0,1}) {
      std::complex<float> value{1,0};
      SpectralParams p;p.primary_mode=PrimaryMode::Table;
      rejects([&]{select_spectral(opt)(&value,nullptr,1,0,p);});
      p={};p.enhancement.sharpen=1;
      rejects([&]{select_spectral(opt)(&value,nullptr,1,0,p);});
      p={};p.enhancement.dehalo=1;
      rejects([&]{select_spectral(opt)(&value,nullptr,1,0,p);});
      const std::complex<float>* slots[]={&value};
      rejects([&]{select_fft3d_temporal(opt)(slots,0,0,1,0,nullptr,{},0,&value);});
      rejects([&]{select_fft3d_temporal(opt)(slots,1,1,1,0,nullptr,{},0,&value);});
      NoisePower noise;noise.mode=PrimaryMode::Table;
      rejects([&]{select_fft3d_temporal(opt)(slots,1,0,1,0,nullptr,noise,0,&value);});
      p={};p.type=4;p.a=1;p.low=p.high=1e20f;
      rejects([&]{select_spectral(opt,true)(&value,nullptr,1,0,p);});
    }
    const auto optimized = select_spectral(0);
    optimized(nullptr, nullptr, 0, 0, {});
    for (int type = -1; type <= 4; ++type)
      for (float exponent : {.5f, 1.0f, 2.0f, .50005f})
        for (std::size_t count : {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 257})
          for (bool mean : {false, true})
          for (bool table : {false, true}) {
            Guarded samples(count), grid(count), primary_storage(count);
            auto* primary = reinterpret_cast<float*>(primary_storage.data) + count;
            for (std::size_t i = 0; i < count; ++i) primary[i] = .25f + float(i % 13);
            auto expected = buffer<std::complex<float>>(count), original_grid = expected;
            for (std::size_t i = 0; i < count; ++i) {
              samples.data[i] = expected[i] = {float(int(i % 19) - 9), float(int(i % 23) - 11)};
              grid.data[i] = original_grid[i] = {float(i % 5) * .25f, float(i % 3) * .5f};
            }
            SpectralParams p{type, 9, .5f, 25, 100, exponent, .25f};
            if (table) {
              p.primary_mode = PrimaryMode::Table;
              p.primary = {primary, count};
            }
            // Enhancement windows share the guarded table boundary, including odd tails.
            if (type == -1 && table) {
              p.enhancement = {.4f, .2f, 16, 400, 2500, {primary, count}, {primary, count}};
            }
            spectral_scalar(expected.data(), mean ? grid.data : nullptr, count, .5f, p);
            optimized(samples.data, mean ? grid.data : nullptr, count, .5f, p);
            CHECK(std::memcmp(grid.data, original_grid.data(), count * sizeof(std::complex<float>)) == 0);
            for (std::size_t i = 0; i < count; ++i) {
              check_near(samples.data[i].real(), expected[i].real(), 2e-6);
              check_near(samples.data[i].imag(), expected[i].imag(), 2e-6);
            }
          }
    // 1. Pure SIMD (exact vector length, count = 32): verify SIMD detection without scalar tail
    {
      SpectralParams overflow_params;
      overflow_params.type = 2;
      overflow_params.a = 3e38f;

      // grid == nullptr: pure SIMD, index 0 overflows (2 * 3e38 = Inf)
      Guarded pure_simd(32);
      pure_simd.data[0] = {2.0f, 0};
      rejects([&] { optimized(pure_simd.data, nullptr, 32, 0, overflow_params); });

      // grid != nullptr: pure SIMD, index 0 overflows (re * gain + mr = Inf)
      Guarded pure_simd_grid(32), grid(32);
      pure_simd_grid.data[0] = {2.0f, 0};
      rejects([&] { optimized(pure_simd_grid.data, grid.data, 32, 0.5f, overflow_params); });
    }

    // 2. SIMD overflow with scalar tail present (count = 33): tail elements are finite & non-overflowing
    {
      SpectralParams overflow_params;
      overflow_params.type = 2;
      overflow_params.a = 3e38f;

      // grid == nullptr: only index 0 overflows, index 32 (tail) has 0.0f (gain * 0 = 0, finite)
      Guarded overflow_test(33);
      overflow_test.data[0] = {2.0f, 0};
      rejects([&] { optimized(overflow_test.data, nullptr, 33, 0, overflow_params); });

      // grid != nullptr: only index 0 overflows in SIMD, index 32 (tail) is finite
      Guarded overflow_test_grid(33), grid(33);
      overflow_test_grid.data[0] = {2.0f, 0};
      rejects([&] { optimized(overflow_test_grid.data, grid.data, 33, 0.5f, overflow_params); });
    }

    // 3. Huge power overflow: only index 0 has huge amplitude, scalar tail is normal and finite
    {
      Guarded huge_no_grid(33);
      huge_no_grid.data[0] = {1e30f, 0};
      for (std::size_t i = 1; i < 33; ++i)
        huge_no_grid.data[i] = {1.0f, 0};
      rejects([&] { optimized(huge_no_grid.data, nullptr, 33, 0, {}); });

      Guarded huge_grid(33), grid(33);
      huge_grid.data[0] = {1e30f, 0};
      for (std::size_t i = 1; i < 33; ++i)
        huge_grid.data[i] = {1.0f, 0};
      rejects([&] { optimized(huge_grid.data, grid.data, 33, 0.5f, {}); });
    }
    // Spatial Kernel Verification (SIMD vs Scalar)
    {
      const auto opt_spatial = select_spatial(0);
      const auto sc_spatial = select_spatial(1);
      CHECK(spatial_target(0) != nullptr);
      CHECK(spatial_target(1) != nullptr);

      for (int count : {1, 2, 3, 4, 7, 8, 9, 12, 15, 16, 17, 24, 31, 32, 33, 48, 64, 65, 128}) {
        // 1. gather_fft3d
        {
          std::vector<float> src(count + 4), wx_a(count);
          for (int i = 0; i < count + 4; ++i) src[i] = float((i * 7 + 3) % 23 - 11) / 12.0f;
          for (int i = 0; i < count; ++i) wx_a[i] = float((i * 5 + 1) % 13) / 13.0f;
          const float wy = 0.75f;
          std::vector<float> blk_opt(count), blk_sc(count);
          opt_spatial.gather_fft3d(src.data() + 1, wx_a.data(), wy, blk_opt.data(), count);
          sc_spatial.gather_fft3d(src.data() + 1, wx_a.data(), wy, blk_sc.data(), count);
          for (int i = 0; i < count; ++i) {
            check_near(blk_opt[i], blk_sc[i]);
          }
        }

        // 2. gather_dfttest
        {
          std::vector<float> src(count + 4), h(count);
          for (int i = 0; i < count + 4; ++i) src[i] = float((i * 11 + 5) % 29 - 14) / 15.0f;
          for (int i = 0; i < count; ++i) h[i] = float((i * 13 + 7) % 19) / 19.0f;
          std::vector<float> blk_opt(count), blk_sc(count);
          opt_spatial.gather_dfttest(src.data() + 2, h.data(), blk_opt.data(), count);
          sc_spatial.gather_dfttest(src.data() + 2, h.data(), blk_sc.data(), count);
          for (int i = 0; i < count; ++i) {
            check_near(blk_opt[i], blk_sc[i]);
          }
        }

        // 3. scatter_fft3d_block
        {
          std::vector<float> inv(count), wx_s(count);
          for (int i = 0; i < count; ++i) inv[i] = float((i * 9 + 4) % 17 - 8) / 9.0f;
          for (int i = 0; i < count; ++i) wx_s[i] = float((i * 3 + 2) % 7) / 7.0f;
          std::vector<float> r_opt(count + 4), r_sc(count + 4);
          for (int i = 0; i < count + 4; ++i) r_opt[i] = r_sc[i] = float((i * 17 - 50) % 23 - 11) / 12.0f;
          opt_spatial.scatter_fft3d_block(inv.data(), wx_s.data(), r_opt.data() + 1, count);
          sc_spatial.scatter_fft3d_block(inv.data(), wx_s.data(), r_sc.data() + 1, count);
          for (int i = 0; i < count + 4; ++i) {
            check_near(r_opt[i], r_sc[i]);
          }
        }

        // 4. scatter_fft3d_row
        {
          std::vector<float> r(count);
          for (int i = 0; i < count; ++i) r[i] = float((i * 13 - 25) % 31 - 15) / 16.0f;
          const float wy = 1.25f;
          std::vector<float> a_opt(count), a_sc(count);
          for (int i = 0; i < count; ++i) a_opt[i] = a_sc[i] = float((i * 7 + 11) % 29 - 14) / 15.0f;
          opt_spatial.scatter_fft3d_row(r.data(), wy, a_opt.data(), count);
          sc_spatial.scatter_fft3d_row(r.data(), wy, a_sc.data(), count);
          for (int i = 0; i < count; ++i) {
            check_near(a_opt[i], a_sc[i]);
          }
        }

        // 5. scatter_dfttest
        {
          std::vector<float> inv(count), h_syn(count);
          for (int i = 0; i < count; ++i) inv[i] = float((i * 19 - 30) % 37 - 18) / 19.0f;
          for (int i = 0; i < count; ++i) h_syn[i] = float((i * 23 + 5) % 41) / 41.0f;
          std::vector<float> a_opt(count + 4), a_sc(count + 4);
          for (int i = 0; i < count + 4; ++i) a_opt[i] = a_sc[i] = float((i * 29 - 40) % 43 - 21) / 22.0f;
          opt_spatial.scatter_dfttest(inv.data(), h_syn.data(), a_opt.data() + 3, count);
          sc_spatial.scatter_dfttest(inv.data(), h_syn.data(), a_sc.data() + 3, count);
          for (int i = 0; i < count + 4; ++i) {
            check_near(a_opt[i], a_sc[i]);
          }
        }

        // 6. store_output_float (FFT3D and DFTTest)
        {
          std::vector<float> accum(count);
          for (int i = 0; i < count; ++i) accum[i] = float((i * 2 - count) % 31 - 15) / 32.0f;
          std::vector<float> dst_opt(count), dst_sc(count);
          opt_spatial.store_output_float(accum.data(), dst_opt.data(), count, true, 1.0f);
          sc_spatial.store_output_float(accum.data(), dst_sc.data(), count, true, 1.0f);
          for (int i = 0; i < count; ++i) {
            check_near(dst_opt[i], dst_sc[i]);
          }
          opt_spatial.store_output_float(accum.data(), dst_opt.data(), count, false, 1.0f / 255.0f);
          sc_spatial.store_output_float(accum.data(), dst_sc.data(), count, false, 1.0f / 255.0f);
          for (int i = 0; i < count; ++i) {
            check_near(dst_opt[i], dst_sc[i]);
          }
        }

        // 7. store_output_uint8 (FFT3D and DFTTest)
        {
          std::vector<float> accum(count);
          for (int i = 0; i < count; ++i) accum[i] = float(i * 5 - 20);
          std::vector<std::uint8_t> dst_opt(count), dst_sc(count);
          opt_spatial.store_output_uint8(accum.data(), dst_opt.data(), count, true, 128.0f, 1.0f, 255.0f);
          sc_spatial.store_output_uint8(accum.data(), dst_sc.data(), count, true, 128.0f, 1.0f, 255.0f);
          CHECK(dst_opt == dst_sc);

          opt_spatial.store_output_uint8(accum.data(), dst_opt.data(), count, false, 0.0f, 1.0f, 255.0f);
          sc_spatial.store_output_uint8(accum.data(), dst_sc.data(), count, false, 0.0f, 1.0f, 255.0f);
          CHECK(dst_opt == dst_sc);
        }

        // 8. store_output_uint16 (FFT3D and DFTTest)
        {
          std::vector<float> accum(count);
          for (int i = 0; i < count; ++i) accum[i] = float(i * 20 - 100);
          std::vector<std::uint16_t> dst_opt(count), dst_sc(count);
          opt_spatial.store_output_uint16(accum.data(), dst_opt.data(), count, true, 512.0f, 4.0f, 1023.0f);
          sc_spatial.store_output_uint16(accum.data(), dst_sc.data(), count, true, 512.0f, 4.0f, 1023.0f);
          CHECK(dst_opt == dst_sc);

          opt_spatial.store_output_uint16(accum.data(), dst_opt.data(), count, false, 0.0f, 4.0f, 1023.0f);
          sc_spatial.store_output_uint16(accum.data(), dst_sc.data(), count, false, 0.0f, 4.0f, 1023.0f);
          CHECK(dst_opt == dst_sc);
        }
      }
    }

    // FFT3D Temporal Kernel Verification (SIMD vs Scalar)
    {
      const auto opt_temporal = select_fft3d_temporal(0);
      const auto sc_temporal = select_fft3d_temporal(1);
      CHECK(fft3d_temporal_target(0) != nullptr);
      CHECK(fft3d_temporal_target(1) != nullptr);

      for (int T : {1, 2, 3, 4, 5}) {
        const int c = T / 2;
        for (int bins : {1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33, 64, 65, 128, 144, 544}) {
          for (bool with_grid : {false, true}) {
            std::vector<std::vector<std::complex<float>>> frames(T, std::vector<std::complex<float>>(bins));
            std::vector<const std::complex<float>*> ptrs(T);
            for (int j = 0; j < T; ++j) {
              for (int k = 0; k < bins; ++k) {
                frames[j][k] = {float((j * 17 + k * 11) % 43 - 21) / 10.0f,
                                float((j * 13 + k * 19) % 37 - 18) / 10.0f};
              }
              ptrs[j] = frames[j].data();
            }

            std::vector<std::complex<float>> grid_data(bins);
            if (with_grid) {
              for (int k = 0; k < bins; ++k) {
                grid_data[k] = {float(k % 7 + 1) * 0.5f, float(k % 5) * 0.25f};
              }
            }

            std::vector<std::complex<float>> out_opt(bins), out_sc(bins);
            const float degrid = with_grid ? 0.75f : 0.0f;
            const float noise = 2.5f;
            const float lower = 0.25f;

            opt_temporal(ptrs.data(), T, c, bins, degrid, with_grid ? grid_data.data() : nullptr, noise, lower, out_opt.data());
            sc_temporal(ptrs.data(), T, c, bins, degrid, with_grid ? grid_data.data() : nullptr, noise, lower, out_sc.data());

            for (int k = 0; k < bins; ++k) {
              if (std::abs(out_opt[k].real() - out_sc[k].real()) > 1e-4f ||
                  std::abs(out_opt[k].imag() - out_sc[k].imag()) > 1e-4f) {
                std::cout << "Mismatch T=" << T << " bins=" << bins << " grid=" << with_grid
                          << " k=" << k << " opt=(" << out_opt[k].real() << "," << out_opt[k].imag()
                          << ") sc=(" << out_sc[k].real() << "," << out_sc[k].imag() << ")\n";
              }
              check_near(out_opt[k].real(), out_sc[k].real(), 1e-4f);
              check_near(out_opt[k].imag(), out_sc[k].imag(), 1e-4f);
            }
          }
        }
      }
    }

    std::cout << "spectral target: " << spectral_target(0) << "; spatial target: " << spatial_target(0)
              << "; temporal target: " << fft3d_temporal_target(0)
              << "; scalar comparison and protected-page tails passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
