#include "algorithms/plan.hpp"
#include "../test.hpp"
#include <future>
using namespace neo_fft;
void windows() {
  auto x = fft3d_axis(64, 8, 2), y = fft3d_axis(48, 6, 2);
  CHECK(x.count == 13 && x.cover == 80 && x.offset == 6);
  CHECK(y.count == 14 && y.cover == 58 && y.offset == 4);
  rejects([] { fft3d_axis(8, 8, 0); });
  auto d = dft_axis(17, 8, 1, 6);
  CHECK(d.cover == 36 && d.offset == 9 && d.count == 15);
  rejects([] { dft_axis(17, 8, 1, 5); });
  CHECK(dft_axis(1, 1, 0, 12).cover == 1);
  rejects([] { dft_axis(1, 1, 1, 0); });
  for (int type = 0; type < 3; ++type)
    for (int b : {4, 8, 9})
      for (int o : {0, 1, b / 2}) {
        const auto g = fft3d_axis(64, b, o);
        const auto win = fft3d_window(b, o, type);
        for (int pos = g.offset; pos < g.offset + g.length; ++pos) {
          double weight = 0;
          for (int i = 0; i < g.count; ++i) {
            const int local = pos - i * g.step;
            if (local >= 0 && local < b)
              weight += double(win.analysis[local]) * win.synthesis[local];
          }
          check_near(weight, 1, 5e-7);
        }
      }
  check_near(dft_raw_window(9, 0, 2, 2.5f), .74, 1e-12);
  check_near(dft_raw_window(9, 1, 2, 2.5f), .50, 1e-12);
  auto w = dft_window(2, 1, 1, 7, 7, 2.5f, 2.5f);
  for (float v : w.h)
    check_near(v, .25, 1e-8);
  check_near(w.wscale, 4, 1e-8);
  for (int id = 0; id <= 11; ++id)
    for (int o : {0, 4, 6}) {
      auto win = dft_window(8, o, 1, id, 7, 2.5f, 2.5f);
      for (int py = 0; py < 8 - o; ++py)
        for (int px = 0; px < 8 - o; ++px) {
          double weight = 0;
          for (int yy = py; yy < 8; yy += 8 - o)
            for (int xx = px; xx < 8; xx += 8 - o)
              weight += double(win.h[yy * 8 + xx]) * win.h[yy * 8 + xx] * 64;
          check_near(weight, 1, 3e-7);
        }
    }
  for (int id = 0; id <= 11; ++id) {
    const auto win = dft_window(3, 0, 0, 7, id, 2.5f, 2.5f);
    const auto tw = dft_raw_window(id, 0, 1, 2.5f);
    check_near(9.0 * win.h[4] * win.h[4], tw * tw, 3e-7);
  }
  rejects([] { dft_raw_window(4, 0, 8, std::numeric_limits<float>::max()); });
  DFTConfig c;
  c.block = 4;
  c.overlap = 0;
  c.swin = 6;
  rejects([&] { Plan p(64, 48, {8, false, false}, c); });

  // 3D DFT window tests
  {
    auto w3d = dft_window_3d(3, 2, 1, 1, 7, 7, 2.5f, 2.5f);
    check_near(w3d.wscale, 4.0f, 1e-6);

    auto w3d_hann = dft_window_3d(3, 2, 1, 1, 7, 0, 2.5f, 2.5f);
    check_near(w3d_hann.wscale, 32.0f / 3.0f, 1e-4);

    auto w1d = dft_window_3d(1, 2, 1, 1, 7, 7, 2.5f, 2.5f);
    check_near(w1d.wscale, 4.0f, 1e-6);

    // Window 6 (flat top) non-unit center gain
    auto w3d_flat = dft_window_3d(3, 2, 1, 1, 7, 6, 2.5f, 2.5f);
    const double tw_c = 0.2810639 + 0.5208972 + 0.1980399; // 1.000001
    check_near(tw_c, 1.000001, 1e-7);
  }
}
void filters() {
  std::complex<float> x{10, 0};
  SpectralParams p;
  p.a = 64;
  spectral_scalar(&x, nullptr, 1, 0, p);
  check_near(x.real(), 3.6);
  x = {10, 0};
  p.floor = .5f;
  spectral_scalar(&x, nullptr, 1, 0, p);
  check_near(x.real(), 5);
  x = {3, 4};
  p = {0, 9, 0, 0, 0, 1, 0};
  spectral_scalar(&x, nullptr, 1, 0, p);
  check_near(x.real(), 1.92);
  check_near(x.imag(), 2.56);
  x = {3, 4};
  p.exponent = .5f;
  spectral_scalar(&x, nullptr, 1, 0, p);
  check_near(x.real(), 2.4);
  check_near(x.imag(), 3.2);
  for (float threshold : {std::nextafter(25.0f, 0.0f), 25.0f, std::nextafter(25.0f, 30.0f)}) {
    x = {3, 4};
    p = {1, threshold};
    spectral_scalar(&x, nullptr, 1, 0, p);
    CHECK(x == (threshold <= 25 ? std::complex<float>(3, 4) : std::complex<float>(0, 0)));
  }
  for (float v : {4.0f, 5.0f, 10.0f, 11.0f}) {
    x = {v, 0};
    p = {3, 2, .5f, 25, 100};
    spectral_scalar(&x, nullptr, 1, 0, p);
    check_near(x.real(), v * (v >= 5 && v <= 10 ? 2 : .5));
  }
  x = {3, 4};
  p = {4, 1, 0, 0, 100};
  spectral_scalar(&x, nullptr, 1, 0, p);
  check_near(x.real(), 3 * std::sqrt(.8));
  for (float exponent : {.49994f, .49995f, .49996f, .50004f, .50005f, .50006f, .99994f, .99995f, .99996f, 1.00004f,
                         1.00005f, 1.00006f, 2.0f}) {
    x = {3, 4};
    p = {0, 9, 0, 0, 0, exponent};
    spectral_scalar(&x, nullptr, 1, 0, p);
    const float a = 16.0f / 25.0f;
    const float gain = std::abs(exponent - 1) < .00005f     ? a
                       : std::abs(exponent - .5f) < .00005f ? std::sqrt(a)
                                                            : std::pow(a, exponent);
    CHECK(x == std::complex<float>(gain * 3, gain * 4));
  }
  x = {std::numeric_limits<float>::max(), 0};
  rejects([&] { spectral_scalar(&x, nullptr, 1, 0, p); });
}
template <class T>
void identity(SampleFormat format) {
  constexpr int w = 35, h = 29, stride = 39;
  auto in = buffer<T>(stride * h), out = buffer<T>(stride * h);
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x) {
      if constexpr (std::is_same_v<T, float>)
        in[y * stride + x] = float((x * 11 + y * 7) % 29) / 20 - .25f;
      else
        in[y * stride + x] = T(((x * 11 + y * 7) % 251) << (format.bits - 8));
    }
  const auto src = checked_plane<const T>(in.data(), w, h, stride * sizeof(T), in.size() * sizeof(T));
  auto dst = checked_plane(out.data(), w, h, stride * sizeof(T), out.size() * sizeof(T));
  for (int type = 0; type < 3; ++type) {
    FFT3DConfig c;
    c.bw = 8;
    c.bh = 6;
    c.ow = 4;
    c.oh = 2;
    c.wintype = type;
    c.sigma = 0;
    Plan plan(w, h, format, c);
    std::fill(out.begin(), out.end(), T(7));
    plan.process(src, dst);
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < stride; ++x) {
        if (x >= w)
          CHECK(out[y * stride + x] == T(7));
        else if constexpr (std::is_same_v<T, float>)
          check_near(out[y * stride + x], std::clamp(in[y * stride + x], 0.0f, 1.0f));
        else
          CHECK(out[y * stride + x] == in[y * stride + x]);
      }
  }
  for (int mode : {0, 1}) {
    DFTConfig c;
    c.block = mode ? 8 : 3;
    c.overlap = mode ? 6 : 12;
    c.mode = mode;
    c.ftype = 2;
    c.sigma = 1;
    c.swin = 7;
    c.zmean = false;
    Plan plan(w, h, format, c);
    plan.process(src, dst);
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        if constexpr (std::is_same_v<T, float>)
          check_near(out[y * stride + x], in[y * stride + x]);
        else
          CHECK(out[y * stride + x] == in[y * stride + x]);
    const auto first = out;
    std::vector<std::future<std::vector<T>>> jobs;
    for (int k = 0; k < 4; ++k)
      jobs.emplace_back(std::async(std::launch::async, [&] {
        auto temp = buffer<T>(stride * h);
        std::fill(temp.begin(), temp.end(), T(7));
        plan.process(src, checked_plane(temp.data(), w, h, stride * sizeof(T), temp.size() * sizeof(T)));
        return temp;
      }));
    for (auto& job : jobs)
      CHECK(job.get() == first);
  }
}

void multithreaded_size_switching() {
  const std::vector<int> sizes = {8, 12, 16, 20, 24, 32, 40, 48, 64};
  std::vector<std::future<void>> futures;
  for (int t = 0; t < 4; ++t) {
    futures.emplace_back(std::async(std::launch::async, [sizes, t] {
      for (int iter = 0; iter < 30; ++iter) {
        for (int b : sizes) {
          FFT3DConfig cfg;
          cfg.bh = b;
          cfg.bw = b;
          cfg.oh = b / 2;
          cfg.ow = b / 2;
          cfg.sigma = 2.0f;
          Plan plan(128, 128, {8, false, false}, cfg);
          auto src_buf = buffer<std::uint8_t>(128 * 128);
          const auto val = std::uint8_t((t + iter + b) % 256);
          std::fill(src_buf.begin(), src_buf.end(), val);
          auto dst_buf = buffer<std::uint8_t>(128 * 128);
          plan.process(checked_plane(src_buf.data(), 128, 128, 128, 128 * 128),
                       checked_plane(dst_buf.data(), 128, 128, 128, 128 * 128));
          CHECK(dst_buf == src_buf);
        }
      }
    }));
  }
  for (auto& f : futures) {
    f.get();
  }
}

int main() {
  try {
    windows();
    filters();
    identity<std::uint8_t>({8, false, false});
    identity<std::uint16_t>({10, false, true});
    identity<std::uint16_t>({16, false, false});
    identity<float>({32, true, true});
    multithreaded_size_switching();
    std::cout << "algorithms: geometry, windows, spectral branches, reconstruction, concurrency passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
