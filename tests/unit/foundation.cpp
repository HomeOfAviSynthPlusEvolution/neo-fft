#include "spectral/fft.hpp"
#include "../reference/dft.hpp"
#include "../test.hpp"
#include <algorithm>
using namespace neo_fft;
int main() {
  try {
    rejects([]{RealFFT3D too_tall(65537,65537,1,FftProfile::scalar);});
    // Churn more than 16 cached lengths, then repeat; eviction must retain
    // shared plans while in use and return valid newly reconstructed plans.
    for(int pass=0;pass<2;++pass) for(int width=17;width<50;++width) {
      RealFFT fft(1,width,FftProfile::scalar);
      std::vector<float> input(width),output(width);input[1]=1;
      std::vector<std::complex<float>> bins(fft.bins());
      fft.forward(input.data(),bins.data());fft.inverse(bins.data(),output.data());
      for(int i=0;i<width;++i)check_near(output[i],input[i],2e-6);
    }
    std::uint16_t a[24]{};
    auto p = checked_plane(a, 5, 3, 16, 42);
    CHECK(p.stride() == 8 && p.stride_bytes() == 16);
    rejects([&] { checked_plane(a, 5, 3, 16, 41); });
    rejects([&] { checked_plane(a, 5, 3, 11, 48); });
    rejects([&] { checked_plane(a, 5, 3, 8, 48); });
    rejects([&] { checked_subplane(p, 4, 0, 2, 1); });
    rejects([&] { mul_size(SIZE_MAX, 2); });
    rejects([&] { add_size(SIZE_MAX, 1); });
    for (int h : {1, 2, 3, 6, 7})
      for (int w : {1, 2, 3, 4, 8, 9}) {
        RealFFT fft(h, w);
        const int k = fft.columns(), rs = w + 3, ss = k + 2;
        const std::size_t rd = rs * h + 5, sd = ss * h + 3;
        for (std::size_t active : {0, 1, 3, 4}) {
          auto in = buffer<float>(rd * 4), out = buffer<float>(rd * 4);
          std::fill(in.begin(), in.end(), -999.0f);
          std::fill(out.begin(), out.end(), -777.0f);
          auto spec = buffer<std::complex<float>>(sd * 4);
          std::fill(spec.begin(), spec.end(), std::complex<float>(-888.0f, 444.0f));
          for (std::size_t b = 0; b < active; ++b)
            for (int y = 0; y < h; ++y)
              for (int x = 0; x < w; ++x)
                in[b * rd + y * rs + x] = float((x * 13 + y * 7 + int(b) * 5) % 19 - 9) / 8;
          const auto source = in;
          BatchLayout r{std::size_t(rs), rd, 4, active}, s{std::size_t(ss), sd, 4, active};
          fft.forward(in.data(), r, spec.data(), s);
          const auto original_spectrum = spec;
          fft.inverse(spec.data(), s, out.data(), r);
          CHECK(in == source && spec == original_spectrum);
          for (std::size_t b = 0; b < 4; ++b) {
            double spatial = 0, spectral = 0;
            const auto oracle = direct_dft(in.data() + b * rd, h, w, rs);
            for (std::size_t i = 0; i < rd; ++i) {
              if (b < active && i / rs < std::size_t(h) && i % rs < std::size_t(w)) {
                check_near(out[b * rd + i], in[b * rd + i]);
                spatial += double(in[b * rd + i]) * in[b * rd + i];
              } else
                CHECK(out[b * rd + i] == -777);
            }
            for (std::size_t i = 0; i < sd; ++i) {
              if (b < active && i / ss < std::size_t(h) && i % ss < std::size_t(k)) {
                auto v = spec[b * sd + i];
                auto d = oracle[(i / ss) * k + i % ss];
                check_near(v.real(), d.real(), 3e-5);
                check_near(v.imag(), d.imag(), 3e-5);
                const int x = int(i % ss);
                spectral += std::norm(std::complex<double>(v)) * (x == 0 || (w % 2 == 0 && x == w / 2) ? 1 : 2);
              } else
                CHECK(spec[b * sd + i] == std::complex<float>(-888, 444));
            }
            if (b < active)
              check_near(spatial, spectral / (w * h), 5e-5);
          }
        }
      }
    RealFFT fft(1, 4);
    fft.forward(nullptr, {0, 0, 4, 0}, nullptr, {0, 0, 4, 0});
    float pulse[]{0, 1, 0, 0}, back[4]{};
    std::complex<float> s[3];
    fft.forward(pulse, s);
    check_near(s[1].imag(), -1);
    fft.inverse(s, back);
    check_near(back[1], 1.0f);

    for (auto profile : {FftProfile::scalar, FftProfile::native
#if NEO_FFT_FFT_X86_TARGETS
                         ,
                         FftProfile::sse2, FftProfile::avx2, FftProfile::avx512
#endif
         }) {
      CHECK(fft_profile_name(profile) != nullptr);
      CHECK(fft_backend_name(profile) != nullptr);
      if (!fft_profile_supported(profile)) {
        rejects([&] { RealFFT(2, 4, profile); });
        continue;
      }
      RealFFT fft_prof(2, 4, profile);
      CHECK(fft_prof.lanes() >= 1);
      float pulse_prof[]{0, 1, 0, 0, 0, 0, 0, 0}, back_prof[8]{};
      std::complex<float> spec_prof[6]{};
      fft_prof.forward(pulse_prof, spec_prof);
      fft_prof.inverse(spec_prof, back_prof);
      check_near(back_prof[1], 1.0f);
    }

    // Dedicated 16x16 Codelet regression tests
    {
      RealFFT fft16(16, 16);
      CHECK(fft16.columns() == 9);

      // 1. 256 unit delta impulses against independent double-precision direct_dft
      for (int py = 0; py < 16; ++py) {
        for (int px = 0; px < 16; ++px) {
          float in[256]{};
          in[py * 16 + px] = 1.0f;
          std::complex<float> spec[16 * 9]{};
          fft16.forward(in, spec);
          const auto oracle = direct_dft(in, 16, 16, 16);
          for (int i = 0; i < 16 * 9; ++i) {
            check_near(spec[i].real(), oracle[i].real(), 5e-6);
            check_near(spec[i].imag(), oracle[i].imag(), 5e-6);
          }
        }
      }

      // Random block test against direct_dft
      float rnd_in[256];
      for (int i = 0; i < 256; ++i) {
        rnd_in[i] = float((i * 37 + 17) % 101 - 50) / 25.0f;
      }
      std::complex<float> rnd_spec[16 * 9]{};
      fft16.forward(rnd_in, rnd_spec);
      const auto oracle_rnd = direct_dft(rnd_in, 16, 16, 16);
      for (int i = 0; i < 16 * 9; ++i) {
        check_near(rnd_spec[i].real(), oracle_rnd[i].real(), 3e-5);
        check_near(rnd_spec[i].imag(), oracle_rnd[i].imag(), 3e-5);
      }

      // 2. Cross inverse transform test (PocketFFT scalar forward -> Codelet inverse)
      {
        RealFFT fft_scalar(16, 16, FftProfile::scalar);
        float orig[256], restored[256]{};
        for (int i = 0; i < 256; ++i) {
          orig[i] = float((i * 43 + 19) % 256 - 128) / 32.0f;
        }
        std::complex<float> scalar_spec[16 * 9]{};
        fft_scalar.forward(orig, scalar_spec);

        fft16.inverse(scalar_spec, restored);
        for (int i = 0; i < 256; ++i) {
          check_near(restored[i], orig[i], 1e-6);
        }
      }

      // 3. Strides, batches (0, 1, 7, 8, 9, 16, 32, 64), input preservation, output sentinels
      for (int rs : {16, 20, 24, 32}) {
        for (int ss : {9, 12, 16}) {
          for (std::size_t active : {0, 1, 7, 8, 9, 16, 32, 64}) {
            const std::size_t total_blocks = 64;
            const std::size_t in_dist = std::size_t(rs) * 16 + 8;
            const std::size_t out_dist = std::size_t(ss) * 16 + 8;

            auto in_buf = buffer<float>(in_dist * total_blocks);
            auto out_buf = buffer<float>(in_dist * total_blocks);
            auto spec_buf = buffer<std::complex<float>>(out_dist * total_blocks);

            std::fill(in_buf.begin(), in_buf.end(), -999.0f);
            std::fill(out_buf.begin(), out_buf.end(), -777.0f);
            std::fill(spec_buf.begin(), spec_buf.end(), std::complex<float>(-888.0f, 444.0f));

            for (std::size_t b = 0; b < active; ++b) {
              for (int y = 0; y < 16; ++y) {
                for (int x = 0; x < 16; ++x) {
                  in_buf[b * in_dist + y * rs + x] = float((x * 17 + y * 13 + int(b) * 7) % 31 - 15) / 16.0f;
                }
              }
            }

            const auto in_orig = in_buf;
            BatchLayout r_layout{std::size_t(rs), in_dist, total_blocks, active};
            BatchLayout s_layout{std::size_t(ss), out_dist, total_blocks, active};

            fft16.forward(in_buf.data(), r_layout, spec_buf.data(), s_layout);
            CHECK(in_buf == in_orig);

            const auto spec_orig = spec_buf;
            fft16.inverse(spec_buf.data(), s_layout, out_buf.data(), r_layout);
            CHECK(spec_buf == spec_orig);

            for (std::size_t b = 0; b < total_blocks; ++b) {
              if (b < active) {
                for (int y = 0; y < 16; ++y) {
                  for (int x = 0; x < 16; ++x) {
                    const auto idx = b * in_dist + y * rs + x;
                    check_near(out_buf[idx], in_buf[idx], 1e-6);
                  }
                }
              } else {
                for (std::size_t i = 0; i < in_dist; ++i) {
                  CHECK(out_buf[b * in_dist + i] == -777.0f);
                }
              }
            }
          }
        }
      }
    }
    {
      RealFFT fft8(8, 8);
      CHECK(fft8.columns() == 5);

      // 1. 64 unit delta impulses against independent double-precision direct_dft
      for (int py = 0; py < 8; ++py) {
        for (int px = 0; px < 8; ++px) {
          float in[64]{};
          in[py * 8 + px] = 1.0f;
          std::complex<float> spec[8 * 5]{};
          fft8.forward(in, spec);
          const auto oracle = direct_dft(in, 8, 8, 8);
          for (int i = 0; i < 8 * 5; ++i) {
            check_near(spec[i].real(), oracle[i].real(), 5e-6);
            check_near(spec[i].imag(), oracle[i].imag(), 5e-6);
          }
        }
      }

      // Random block test against direct_dft
      float rnd_in[64];
      for (int i = 0; i < 64; ++i) {
        rnd_in[i] = float((i * 31 + 17) % 97 - 48) / 25.0f;
      }
      std::complex<float> rnd_spec[8 * 5]{};
      fft8.forward(rnd_in, rnd_spec);
      const auto oracle_rnd = direct_dft(rnd_in, 8, 8, 8);
      for (int i = 0; i < 8 * 5; ++i) {
        check_near(rnd_spec[i].real(), oracle_rnd[i].real(), 1e-5);
        check_near(rnd_spec[i].imag(), oracle_rnd[i].imag(), 1e-5);
      }

      // 2. Cross inverse transform test (PocketFFT scalar forward -> Codelet inverse)
      {
        RealFFT fft_scalar(8, 8, FftProfile::scalar);
        float orig[64], restored[64]{};
        for (int i = 0; i < 64; ++i) {
          orig[i] = float((i * 37 + 19) % 256 - 128) / 32.0f;
        }
        std::complex<float> scalar_spec[8 * 5]{};
        fft_scalar.forward(orig, scalar_spec);

        fft8.inverse(scalar_spec, restored);
        for (int i = 0; i < 64; ++i) {
          check_near(restored[i], orig[i], 1e-6);
        }
      }

      for (int rs : {8, 12, 16, 24}) {
        for (int ss : {5, 8, 12}) {
          for (std::size_t active : {0, 1, 7, 8, 9, 16, 32, 64}) {
            const std::size_t total_blocks = 64;
            const std::size_t in_dist = std::size_t(rs) * 8 + 8;
            const std::size_t out_dist = std::size_t(ss) * 8 + 8;

            auto in_buf = buffer<float>(in_dist * total_blocks);
            auto out_buf = buffer<float>(in_dist * total_blocks);
            auto spec_buf = buffer<std::complex<float>>(out_dist * total_blocks);

            std::fill(in_buf.begin(), in_buf.end(), -999.0f);
            std::fill(out_buf.begin(), out_buf.end(), -777.0f);
            std::fill(spec_buf.begin(), spec_buf.end(), std::complex<float>(-888.0f, 444.0f));

            for (std::size_t b = 0; b < active; ++b) {
              for (int y = 0; y < 8; ++y) {
                for (int x = 0; x < 8; ++x) {
                  in_buf[b * in_dist + y * rs + x] = float((x * 19 + y * 11 + int(b) * 5) % 29 - 14) / 15.0f;
                }
              }
            }

            const auto in_orig = in_buf;
            BatchLayout r_layout{std::size_t(rs), in_dist, total_blocks, active};
            BatchLayout s_layout{std::size_t(ss), out_dist, total_blocks, active};

            fft8.forward(in_buf.data(), r_layout, spec_buf.data(), s_layout);
            CHECK(in_buf == in_orig);

            const auto spec_orig = spec_buf;
            fft8.inverse(spec_buf.data(), s_layout, out_buf.data(), r_layout);
            CHECK(spec_buf == spec_orig);

            for (std::size_t b = 0; b < total_blocks; ++b) {
              if (b < active) {
                for (int y = 0; y < 8; ++y) {
                  for (int x = 0; x < 8; ++x) {
                    const auto idx = b * in_dist + y * rs + x;
                    check_near(out_buf[idx], in_buf[idx], 1e-5);
                  }
                }
              } else {
                for (std::size_t i = 0; i < in_dist; ++i) {
                  CHECK(out_buf[b * in_dist + i] == -777.0f);
                }
              }
            }
          }
        }
      }
    }
    {
      RealFFT fft32(32, 32);
      CHECK(fft32.columns() == 17);

      // 1. Impulses across diagonals and boundaries against double-precision direct_dft
      for (int p = 0; p < 32; ++p) {
        float in[1024]{};
        in[p * 32 + p] = 1.0f;
        std::complex<float> spec[32 * 17]{};
        fft32.forward(in, spec);
        const auto oracle = direct_dft(in, 32, 32, 32);
        for (int i = 0; i < 32 * 17; ++i) {
          check_near(spec[i].real(), oracle[i].real(), 1.5e-4);
          check_near(spec[i].imag(), oracle[i].imag(), 1.5e-4);
        }
      }

      // Random block test against direct_dft
      float rnd_in[1024];
      for (int i = 0; i < 1024; ++i) {
        rnd_in[i] = float((i * 47 + 23) % 103 - 51) / 25.0f;
      }
      std::complex<float> rnd_spec[32 * 17]{};
      fft32.forward(rnd_in, rnd_spec);
      const auto oracle_rnd = direct_dft(rnd_in, 32, 32, 32);
      for (int i = 0; i < 32 * 17; ++i) {
        check_near(rnd_spec[i].real(), oracle_rnd[i].real(), 1.5e-4);
        check_near(rnd_spec[i].imag(), oracle_rnd[i].imag(), 1.5e-4);
      }

      // 2. Cross inverse transform test (PocketFFT scalar forward -> Codelet inverse)
      {
        RealFFT fft_scalar(32, 32, FftProfile::scalar);
        float orig[1024], restored[1024]{};
        for (int i = 0; i < 1024; ++i) {
          orig[i] = float((i * 53 + 29) % 256 - 128) / 32.0f;
        }
        std::complex<float> scalar_spec[32 * 17]{};
        fft_scalar.forward(orig, scalar_spec);

        fft32.inverse(scalar_spec, restored);
        for (int i = 0; i < 1024; ++i) {
          check_near(restored[i], orig[i], 1e-5);
        }
      }

      // 3. Strides, batches (0, 1, 7, 8, 16, 32), input preservation, output sentinels
      for (int rs : {32, 40, 48}) {
        for (int ss : {17, 20, 24}) {
          for (std::size_t active : {0, 1, 7, 8, 16, 32}) {
            const std::size_t total_blocks = 32;
            const std::size_t in_dist = std::size_t(rs) * 32 + 8;
            const std::size_t out_dist = std::size_t(ss) * 32 + 8;

            auto in_buf = buffer<float>(in_dist * total_blocks);
            auto out_buf = buffer<float>(in_dist * total_blocks);
            auto spec_buf = buffer<std::complex<float>>(out_dist * total_blocks);

            std::fill(in_buf.begin(), in_buf.end(), -999.0f);
            std::fill(out_buf.begin(), out_buf.end(), -777.0f);
            std::fill(spec_buf.begin(), spec_buf.end(), std::complex<float>(-888.0f, 444.0f));

            for (std::size_t b = 0; b < active; ++b) {
              for (int y = 0; y < 32; ++y) {
                for (int x = 0; x < 32; ++x) {
                  in_buf[b * in_dist + y * rs + x] = float((x * 23 + y * 17 + int(b) * 9) % 37 - 18) / 19.0f;
                }
              }
            }

            const auto in_orig = in_buf;
            BatchLayout r_layout{std::size_t(rs), in_dist, total_blocks, active};
            BatchLayout s_layout{std::size_t(ss), out_dist, total_blocks, active};

            fft32.forward(in_buf.data(), r_layout, spec_buf.data(), s_layout);
            CHECK(in_buf == in_orig);

            const auto spec_orig = spec_buf;
            fft32.inverse(spec_buf.data(), s_layout, out_buf.data(), r_layout);
            CHECK(spec_buf == spec_orig);

            for (std::size_t b = 0; b < total_blocks; ++b) {
              if (b < active) {
                for (int y = 0; y < 32; ++y) {
                  for (int x = 0; x < 32; ++x) {
                    const auto idx = b * in_dist + y * rs + x;
                    check_near(out_buf[idx], in_buf[idx], 1e-5);
                  }
                }
              } else {
                for (std::size_t i = 0; i < in_dist; ++i) {
                  CHECK(out_buf[b * in_dist + i] == -777.0f);
                }
              }
            }
          }
        }
      }
    }

    // 3D Real FFT test
    for (int d : {1, 3}) {
      for (int h : {2, 3}) {
        for (int w : {2, 4, 5}) {
          RealFFT3D fft3d(d, h, w);
          CHECK(fft3d.depth() == d);
          CHECK(fft3d.height() == h);
          CHECK(fft3d.width() == w);
          CHECK(fft3d.columns() == w / 2 + 1);
          CHECK(fft3d.samples() == std::size_t(d) * h * w);
          CHECK(fft3d.bins() == std::size_t(d) * h * (w / 2 + 1));

          const std::size_t n_samples = fft3d.samples();
          const std::size_t n_bins = fft3d.bins();
          std::vector<float> in(n_samples), out(n_samples, -999.0f);
          std::vector<std::complex<float>> spec(n_bins);

          for (std::size_t i = 0; i < n_samples; ++i) {
            in[i] = float(int(i * 17 + 3) % 23 - 11) / 8.0f;
          }

          fft3d.forward(in.data(), spec.data());
          fft3d.inverse(spec.data(), out.data());

          const auto oracle = direct_dft_3d(in.data(), d, h, w);
          double spatial_energy = 0, spectral_energy = 0;
          const int k = w / 2 + 1;
          for (std::size_t i = 0; i < n_samples; ++i) {
            spatial_energy += double(in[i]) * in[i];
          }
          for (int kz = 0; kz < d; ++kz) {
            for (int ky = 0; ky < h; ++ky) {
              for (int kx = 0; kx < k; ++kx) {
                const auto idx = (kz * h + ky) * k + kx;
                check_near(spec[idx].real(), oracle[idx].real(), 5e-5);
                check_near(spec[idx].imag(), oracle[idx].imag(), 5e-5);
                const double weight = (kx == 0 || (w % 2 == 0 && kx == w / 2)) ? 1.0 : 2.0;
                spectral_energy += std::norm(std::complex<double>(spec[idx])) * weight;
              }
            }
          }
          check_near(spatial_energy, spectral_energy / double(n_samples), 5e-5);
          for (std::size_t i = 0; i < n_samples; ++i) {
            check_near(out[i], in[i], 5e-5);
          }
        }
      }
    }

    // Grouped 3D calls: independent forward oracle, partial groups and stride guards.
    for(int T:{3,5})for(int W:{4,5})for(int count:{0,1,7,8,9}) {
      RealFFT3D transform(T,3,W);
      const auto samples=transform.samples(),bins=transform.bins(),rd=samples+5,sd=bins+3;
      std::vector<float> input(rd*10,-91),back(rd*10,-92);
      const std::complex<float> sentinel{-93,94};
      std::vector<std::complex<float>> spectrum(sd*10,sentinel);
      for(int b=0;b<count;++b)for(std::size_t i=0;i<samples;++i)
        input[b*rd+i]=float(int((i*17+b*7)%23)-11)/8.f;
      const auto saved=input;
      transform.forward(input.data(),count,rd,spectrum.data(),sd);
      transform.inverse(spectrum.data(),count,sd,back.data(),rd);
      CHECK(input==saved);
      for(int b=0;b<count;++b) {
        const auto oracle=direct_dft_3d(input.data()+b*rd,T,3,W);
        for(std::size_t i=0;i<bins;++i) {
          check_near(spectrum[b*sd+i].real(),oracle[i].real(),5e-5);
          check_near(spectrum[b*sd+i].imag(),oracle[i].imag(),5e-5);
        }
        for(std::size_t i=0;i<samples;++i)check_near(back[b*rd+i],input[b*rd+i],5e-5);
      }
      for(int b=0;b<10;++b) {
        for(std::size_t i=b<count ? samples : 0;i<rd;++i)CHECK(back[b*rd+i]==-92);
        for(std::size_t i=b<count ? bins : 0;i<sd;++i)CHECK(spectrum[b*sd+i]==sentinel);
      }
    }

    // Common DFT volumes: native packed axes/codelets versus a binary64 DFT.
    for(int T:{3,5})for(int S:{12,16}) {
      const std::size_t samples=T*S*S,bins=T*S*(S/2+1),rd=samples+5,sd=bins+3;
      std::vector<float> input(rd*10,-91);
      for(int b=0;b<10;++b)for(std::size_t i=0;i<samples;++i)
        input[b*rd+i]=float(int((i*17+(b%2)*7)%23)-11)/8.f;
      const auto saved=input;
      const auto oracle0=direct_dft_3d(input.data(),T,S,S),oracle1=direct_dft_3d(input.data()+rd,T,S,S);
      for(auto profile:{FftProfile::scalar,FftProfile::sse2,FftProfile::avx2,FftProfile::avx512,FftProfile::native}) {
        if(!fft_profile_supported(profile))continue;
        RealFFT3D transform(T,S,S,profile);
        CHECK(!transform.try_inverse_center(nullptr,0,0,nullptr,0));
        CHECK(!transform.try_inverse_center(nullptr,1,0,nullptr,0));
        for(int count:{0,1,2,7,8,9}) {
          std::vector<float> back(rd*10,-92);
          const std::complex<float> sentinel{-93,94};
          std::vector<std::complex<float>> spectrum(sd*10,sentinel);
          transform.forward(input.data(),count,rd,spectrum.data(),sd);
          const auto original_spectrum=spectrum;
          transform.inverse(spectrum.data(),count,sd,back.data(),rd);
          CHECK(input==saved && spectrum==original_spectrum);
          for(int b=0;b<count;++b) {
            const auto& oracle=b%2 ? oracle1 : oracle0;
            for(std::size_t i=0;i<bins;++i) {
              check_near(spectrum[b*sd+i].real(),oracle[i].real(),5e-4);
              check_near(spectrum[b*sd+i].imag(),oracle[i].imag(),5e-4);
            }
            for(std::size_t i=0;i<samples;++i)check_near(back[b*rd+i],input[b*rd+i],2e-6);
          }
          for(int b=0;b<10;++b) {
            for(std::size_t i=b<count ? samples : 0;i<rd;++i)CHECK(back[b*rd+i]==-92);
            for(std::size_t i=b<count ? bins : 0;i<sd;++i)CHECK(spectrum[b*sd+i]==sentinel);
          }
          // Independent binary64 forward spectra test selective inverse sign,
          // center index and scale; non-target slices and inactive volumes are guards.
          for(int b=0;b<count;++b) {
            const auto& oracle=b%2 ? oracle1 : oracle0;
            for(std::size_t i=0;i<bins;++i)spectrum[b*sd+i]=std::complex<float>(oracle[i]);
          }
          const auto saved_spectrum=spectrum;
          const auto center=std::size_t(T/2)*S*S;
          std::vector<float> selected(rd*10+2,-96);
          const bool used=transform.try_inverse_center(spectrum.data(),count,sd,selected.data()+1+center,rd);
          CHECK(spectrum==saved_spectrum);
#if defined(NEO_FFT_HAS_AVX2_CODELET)
          if(S==16 && count>1 && transform.lanes()>1 &&
             (profile==FftProfile::avx2 || profile==FftProfile::avx512))CHECK(used);
#endif
          if(count<=1 || S!=16 || profile==FftProfile::scalar || profile==FftProfile::sse2)CHECK(!used);
          CHECK(selected.front()==-96 && selected.back()==-96);
          for(int b=0;b<10;++b)for(std::size_t i=0;i<rd;++i) {
            if(used && b<count && i>=center && i<center+S*S)
              check_near(selected[1+b*rd+i],input[b*rd+i],2e-6);
            else CHECK(selected[1+b*rd+i]==-96);
          }
        }
      }
    }

    // 3D impulse & constant checks
    {
      const int T = 3, H = 2, W = 4;
      RealFFT3D fft3d(T, H, W);
      const auto n_samples = fft3d.samples();
      const auto n_bins = fft3d.bins();
      const int K = fft3d.columns();

      // Unit impulse at z=1, y=0, x=0
      std::vector<float> impulse(n_samples, 0.0f);
      impulse[(1 * H + 0) * W + 0] = 1.0f;
      std::vector<std::complex<float>> spec(n_bins);
      fft3d.forward(impulse.data(), spec.data());

      const double pi = std::acos(-1.0);
      for (int kt = 0; kt < T; ++kt) {
        const double expected_angle = -2 * pi * double(kt) / double(T);
        for (int ky = 0; ky < H; ++ky) {
          for (int kx = 0; kx < K; ++kx) {
            const auto idx = (kt * H + ky) * K + kx;
            check_near(spec[idx].real(), float(std::cos(expected_angle)), 5e-5);
            check_near(spec[idx].imag(), float(std::sin(expected_angle)), 5e-5);
          }
        }
      }

      // Constant q
      const float q = 3.5f;
      std::vector<float> constant_in(n_samples, q), back(n_samples, 0.0f);
      fft3d.forward(constant_in.data(), spec.data());
      check_near(spec[0].real(), float(n_samples) * q, 5e-5);
      check_near(spec[0].imag(), 0.0f, 5e-5);
      fft3d.inverse(spec.data(), back.data());
      for (float v : back) {
        check_near(v, q, 5e-5);
      }
    }

    std::cout << "foundation: checked views, direct DFT, Parseval, strided partial batches, 8x8, 16x16, 32x32 codelets, "
                 "multi-target PocketFFT profiles passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
