#include "spectral/fft.hpp"
#include "../reference/dft.hpp"
#include "../test.hpp"
#include <algorithm>
using namespace neo_fft;
int main() {
  try {
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
          std::fill(in.begin(), in.end(), -999);
          std::fill(out.begin(), out.end(), -777);
          auto spec = buffer<std::complex<float>>(sd * 4);
          std::fill(spec.begin(), spec.end(), std::complex<float>(-888, 444));
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
    s[0] = {1, 1};
    rejects([&] { fft.inverse(s, back); });
    pulse[0] = NAN;
    rejects([&] { fft.forward(pulse, s); });

    for (auto profile : {FftProfile::scalar, FftProfile::native
#if NEO_FFT_FFT_X86_TARGETS
                         ,
                         FftProfile::sse2, FftProfile::avx2, FftProfile::avx512
#endif
         }) {
      RealFFT fft_prof(2, 4, profile);
      CHECK(fft_prof.lanes() >= 1);
      CHECK(fft_profile_name(profile) != nullptr);
      CHECK(fft_backend_name(profile) != nullptr);
      float pulse_prof[]{0, 1, 0, 0, 0, 0, 0, 0}, back_prof[8]{};
      std::complex<float> spec_prof[6]{};
      fft_prof.forward(pulse_prof, spec_prof);
      fft_prof.inverse(spec_prof, back_prof);
      check_near(back_prof[1], 1.0f);
    }

    std::cout << "foundation: checked views, direct DFT, Parseval, strided partial batches, multi-target PocketFFT "
                 "profiles passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
