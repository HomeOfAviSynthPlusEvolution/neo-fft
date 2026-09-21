#include "algorithms/windows.hpp"
#include <array>

namespace neo_fft {
AxisWindow fft3d_window(int block, int overlap, int type) {
  require(block >= 2 && overlap >= 0 && overlap <= block / 2, "FFT3D window geometry invalid");
  require(type >= 0 && type <= 2, "FFT3D wintype outside 0..2");
  AxisWindow win{buffer<float>(block), buffer<float>(block)};
  std::fill(win.analysis.begin(), win.analysis.end(), 1.0f);
  std::fill(win.synthesis.begin(), win.synthesis.end(), 1.0f);
  constexpr float pi = 3.14159265358979323846f;
  for (int j = 0; j < overlap; ++j) {
    const float left = std::cos(pi * (float(j - overlap) + 0.5f) / (2.0f * float(overlap)));
    const float right = std::cos(pi * (float(j) + 0.5f) / (2.0f * float(overlap)));
    for (int side = 0; side < 2; ++side) {
      const float c = side ? right : left;
      const int pos = side ? block - overlap + j : j;
      const float a = type == 0 ? c : type == 1 ? std::sqrt(c) : 1.0f;
      win.analysis[pos] = a;
      win.synthesis[pos] = type == 0 ? c : type == 1 ? (a * a) * a : c * c;
    }
  }
  return win;
}
namespace {
double i0a(double arg) {
  double num = 1, den = 1, total = 1;
  for (int k = 1; k <= 14; ++k) {
    num *= arg / 2;
    den *= k;
    const double v = num / den;
    total += v * v;
    require(std::isfinite(total), "DFTTest Kaiser window overflow");
    if (v <= 1e-8)
      break;
  }
  return total;
}
} // namespace
double dft_raw_window(int id, int j, int length, float beta) {
  require(id >= 0 && id <= 11 && length > 0 && j >= 0 && j < length && std::isfinite(beta) && beta >= 0,
          "DFTTest window parameter invalid");
  constexpr double pi = 3.1415926535897932384626433832795;
  const double n = j + 0.5, r = n / length, theta = 2 * pi * r;
  auto c = [&](int m) {
    return std::cos(m * theta);
  };
  double value = 0;
  switch (id) {
    case 0:
      value = 0.5 - 0.5 * c(1);
      break;
    case 1:
      value = 0.53836 - 0.46164 * c(1);
      break;
    case 2:
      value = 0.42 - 0.5 * c(1) + 0.08 * c(2);
      break;
    case 3:
      value = 0.35875 - 0.48829 * c(1) + 0.14128 * c(2) - 0.01168 * c(3);
      break;
    case 4:
      value = i0a(pi * beta * std::sqrt(1 - (2 * r - 1) * (2 * r - 1))) / i0a(pi * beta);
      break;
    case 5: {
      constexpr std::array<double, 7> a{0.27105140069342415,   -0.433297939234486060, 0.218122999543110620,
                                        -0.065925446388030898, 0.010811742098372268,  -7.7658482522509342e-4,
                                        1.3887217350903198e-5};
      value = a[0];
      for (int m = 1; m < 7; ++m)
        value += a[m] * c(m);
      break;
    }
    case 6:
      value = 0.2810639 - 0.5208972 * c(1) + 0.1980399 * c(2);
      break;
    case 7:
      value = 1;
      break;
    case 8:
      value = (2.0 / length) * (length / 2.0 - std::abs(n - length / 2.0));
      break;
    case 9:
      value = 0.62 - 0.48 * (r - 0.5) - 0.38 * c(1);
      break;
    case 10:
      value = 0.355768 - 0.487396 * c(1) + 0.144232 * c(2) - 0.012604 * c(3);
      break;
    case 11:
      value = 0.3635819 - 0.4891775 * c(1) + 0.1365995 * c(2) - 0.0106411 * c(3);
      break;
  }
  require(std::isfinite(value), "DFTTest non-finite window");
  return value;
}
DftWindow dft_window_3d(int tbsize, int block, int overlap, int mode, int spatial, int temporal, float sbeta, float tbeta) {
  dimension(tbsize);
  dimension(block);
  require(tbsize > 0 && tbsize <= 15 && tbsize % 2 == 1, "DFTTest tbsize must be odd integer in 1..15");
  require(mode == 0 || mode == 1, "DFTTest window mode invalid");
  require(overlap >= 0 && overlap < block, "DFTTest window overlap invalid");
  auto raw = buffer<double>(block), sw = buffer<double>(block);
  for (int j = 0; j < block; ++j)
    raw[j] = dft_raw_window(spatial, j, block, sbeta);
  sw = raw;
  if (mode == 1) {
    const int step = block - overlap;
    for (int j = 0; j < block; ++j) {
      double d = 0;
      for (int k = j; k >= 0; k -= step)
        d += raw[k] * raw[k];
      for (std::int64_t k = std::int64_t(j) + step; k < block; k += step)
        d += raw[k] * raw[k];
      require(std::isfinite(d) && d > 0, "DFTTest zero or invalid window denominator");
      sw[j] = raw[j] / std::sqrt(d);
    }
  }
  auto tw = buffer<double>(tbsize);
  for (int z = 0; z < tbsize; ++z)
    tw[z] = dft_raw_window(temporal, z, tbsize, tbeta);

  const double volume = double(tbsize) * double(block) * double(block);
  const double inv_sqrt_v = 1.0 / std::sqrt(volume);
  const std::size_t total_samples = mul_size(std::size_t(tbsize), mul_size(std::size_t(block), std::size_t(block)));
  DftWindow w{buffer<float>(total_samples), 0};
  float energy = 0.0f;
  std::size_t idx = 0;
  for (int z = 0; z < tbsize; ++z)
    for (int y = 0; y < block; ++y)
      for (int x = 0; x < block; ++x) {
        const double v = tw[z] * sw[y] * sw[x] * inv_sqrt_v;
        require(std::isfinite(v) && std::abs(v) <= std::numeric_limits<float>::max(), "DFTTest window overflow");
        const float h = float(v);
        w.h[idx++] = h;
        energy += h * h;
      }
  require(std::isfinite(energy) && energy > 0, "DFTTest zero or invalid window energy");
  w.wscale = 1.0f / energy;
  require(std::isfinite(w.wscale) && w.wscale > 0, "DFTTest invalid wscale");
  return w;
}

DftWindow dft_window(int block, int overlap, int mode, int spatial, int temporal, float sbeta, float tbeta) {
  return dft_window_3d(1, block, overlap, mode, spatial, temporal, sbeta, tbeta);
}
} // namespace neo_fft
