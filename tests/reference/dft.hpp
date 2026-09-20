#pragma once
#include <complex>
#include <vector>
// Independent binary64 definition. No reference filter or FFT-library calls.
inline std::vector<std::complex<double>> direct_dft(const float* p, int h, int w, int stride) {
  std::vector<std::complex<double>> out(h * (w / 2 + 1));
  const double pi = std::acos(-1.0);
  for (int ky = 0; ky < h; ++ky)
    for (int kx = 0; kx <= w / 2; ++kx)
      for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
          const double angle = -2 * pi * (double(ky) * y / h + double(kx) * x / w);
          out[ky * (w / 2 + 1) + kx] += double(p[y * stride + x]) * std::complex<double>(std::cos(angle), std::sin(angle));
        }
  return out;
}
