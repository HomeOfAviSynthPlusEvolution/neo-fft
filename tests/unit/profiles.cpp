#include "algorithms/plan.hpp"
#include "../test.hpp"
#include <algorithm>

using namespace neo_fft;
namespace {
void curves() {
  DFTCurves c;
  c.shared = {1, 16, 0, 0}; // unsorted, root before interpolation
  auto p = dft_profile(c, 1, 8, 99, 1);
  check_near(p[2 * 5 + 2], 4, 1e-6);
  for (int system : {0, 1})
    for (int T : {1, 3, 5, 15})
      for (int S : {1, 5, 8}) {
        c.system = system;
        p = dft_profile(c, T, S, 99, .25f);
        const int d = (T > 1) + 2 * (S > 1);
        for (int t = 0; t < T; ++t)
          for (int y = 0; y < S; ++y)
            for (int x = 0; x <= S / 2; ++x) {
              const double ft = T == 1 ? 0 : double(std::min(t, T - t)) / (T / 2);
              const double fy = S == 1 ? 0 : double(std::min(y, S - y)) / (S / 2);
              const double fx = S == 1 ? 0 : double(x) / (S / 2);
              double expected = 0;
              if (d) {
                if (system)
                  expected = 16 * std::sqrt((ft * ft + fy * fy + fx * fx) / d);
                else {
                  const double root = std::pow(16, 1.0 / d);
                  expected = (T == 1 ? 1 : root * ft) * (S == 1 ? 1 : root * fy) * (S == 1 ? 1 : root * fx);
                }
              }
              check_near(p[(t * S + y) * (S / 2 + 1) + x], expected * 4, 2e-5);
            }
      }
  c = {};
  c.x = {0, 0, 1, 100};
  c.system = 1;
  p = dft_profile(c, 1, 8, 9, 1);
  for (auto v : p)
    check_near(v, 3); // radial axis-only uses default temporal curve
  c.time = {0, 4, .5f, 16, 1, 36};
  p = dft_profile(c, 1, 8, 9, 1);
  check_near(p[0], 2);
  check_near(p[4 * 5 + 4], 6);
  c.shared = {0, 25, 1, 25};
  p = dft_profile(c, 1, 8, 9, 1);
  for (auto v : p)
    check_near(v, 25);
  c = {};
  c.x = {0, 99, 1, 99};
  check_near(dft_profile(c, 1, 1, 7, 2)[0], 3.5);
  c.time = {0, 13, 1, 17};
  check_near(dft_profile(c, 1, 1, 7, 2)[0], 6.5);
  c.shared = {0, 23, 1, 29};
  check_near(dft_profile(c, 1, 1, 7, 2)[0], 11.5);
  // Exact interior knot, neighboring bins, and arbitrary order.
  c.shared = {1, 16, .5f, 4, 0, 0};
  c.system = 1;
  p = dft_profile(c, 5, 1, 0, 1);
  check_near(p[0], 0);
  check_near(p[1], 4);
  check_near(p[2], 16);
  check_near(p[3], 16);
  for (const auto& invalid : std::vector<std::vector<float>>{
           {0}, {0, 1}, {0, 1, .5f, 1}, {0, 1, 0, 2, 1, 3}, {0, -1, 1, 1}, {-.1f, 1, 1, 1}, {0, 1, 1, INFINITY}}) {
    c.x = invalid;
    rejects([&] { validate(c); }); // ignored axis still validated
  }
  c = {};
  c.system = 2;
  rejects([&] { validate(c); });
  c = {};
  c.shared = {0, 3e38f, 1, 3e38f};
  rejects([&] { dft_profile(c, 1, 1, 0, .1f); });
}
void integration() {
  // Constant curves reproduce scalar sigma for all five filters and both dispatches.
  for (int opt : {0, 1})
    for (int type = 0; type <= 4; ++type)
      for (int T : {1, 3}) {
        DFTConfig c;
        c.block = 5;
        c.overlap = 2;
        c.tbsize = T;
        c.ftype = type;
        c.opt = opt;
        c.sigma = type < 2 ? 4 : .25f;
        Plan uniform(24, 20, {32, true, false}, c);
        c.curves.shared = {0, c.sigma, 1, c.sigma};
        Plan curved(24, 20, {32, true, false}, c);
        std::vector<float> src(24 * 20), a(src.size()), b(src.size());
        for (std::size_t i = 0; i < src.size(); ++i)
          src[i] = float(i % 37) / 40;
        std::vector<span2d::Plane<const float>> slots(T, {src.data(), 24, 20, 24 * sizeof(float)});
        span2d::Span<const span2d::Plane<const float>> input{slots.data(), slots.size()};
        uniform.process(input, {a.data(), 24, 20, 24 * sizeof(float)});
        curved.process(input, {b.data(), 24, 20, 24 * sizeof(float)});
        for (std::size_t i = 0; i < a.size(); ++i)
          check_near(a[i], b[i], 3e-6);
      }
}
} // namespace
int main() {
  try {
    curves();
    integration();
    std::cout << "profiles: independent curve oracles and integration passed\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
