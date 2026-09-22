#include "algorithms/plan.hpp"
#include "../test.hpp"
#include "../reference/dft.hpp"
#include <future>
#include <algorithm>

using namespace neo_fft;
int main() {
  try {
    for (int W : {5, 8, 10}) for (int H : {5, 8, 9}) {
      const std::array<float, 4> sigma{2, 4, 6, 8};
      const auto p = fft3d_profile(W, H, sigma);
      for (int y = 0; y < H; ++y) for (int x = 0; x <= W / 2; ++x) {
        const double fy = double(H - 2 * std::abs(y - H / 2)) / H, fx = double(x) / (W / 2 + 1);
        const double f = std::sqrt((fx * fx + fy * fy) / 2), a = std::sqrt(.5) / 4, b = 2 * a;
        const double s = f < a ? 8 - 2 * f / a : f < b ? 6 - 2 * (f-a) / (b-a) : 2 + 2 * (1-f) / (1-b);
        check_near(p[y * (W / 2 + 1) + x], s * s * W * H, .001);
      }
      EnhancementConfig c;
      c.sharpen = .4f; c.dehalo = .2f;
      const auto e = enhancement_tables(W, H, 4, c);
      auto params = SpectralParams{}; params.type = -2; params.enhancement = e.view(c);
      auto in = buffer<std::complex<float>>(p.size()), expected = in, actual = in;
      double max_halo = 0;
      for (int y = 0; y < H; ++y) for (int x = 0; x <= W / 2; ++x) {
        const double dy = y < H / 2 ? y : H-y;
        const double d2 = dy*dy / (H/2 * (H/2)) + double(x*x) / (W/2 * (W/2));
        max_halo = std::max(max_halo, std::exp(-.7*d2*4) - std::exp(-d2*4));
      }
      for (int y = 0; y < H; ++y) for (int x = 0; x <= W / 2; ++x) {
        const auto k = y*(W/2+1)+x;
        const double dy = y < H/2 ? y : H-y;
        const double d2 = dy*dy / (H/2 * (H/2)) + double(x*x) / (W/2 * (W/2));
        check_near(e.sharpen[k], 1-std::exp(-d2/(2*.3*.3)), 2e-6);
        check_near(e.halo[k], (std::exp(-.7*d2*4)-std::exp(-d2*4))/max_halo, 2e-6);
        in[k] = {float(k+1)*3, float(k%7)*2};
        const double q = std::norm(std::complex<double>(in[k])) + 1e-15;
        const double A = 16.*16*W*H, B = 80.*80*W*H, C = 50.*50*W*H;
        const double S = 1 + c.sharpen*e.sharpen[k]*std::sqrt(q*B/((q+A)*(q+B)));
        const double D = (q+C)/(q+C+c.dehalo*e.halo[k]*q);
        expected[k] = in[k]*float(S*D);
      }
      for (int opt : {0,1}) {
        actual = in;
        select_spectral(opt)(actual.data(), nullptr, actual.size(), 0, params);
        for (std::size_t k = 0; k < actual.size(); ++k) {
          check_near(actual[k].real(), expected[k].real(), 5e-5);
          check_near(actual[k].imag(), expected[k].imag(), 5e-5);
        }
      }
    }
    // Temporal table broadcast uses T*P, independently evaluated with a direct DFT.
    for (int T = 1; T <= 5; ++T) for (int opt : {0,1}) {
      constexpr int bins = 17;
      std::complex<float> inputs[5][bins], output[bins];
      const std::complex<float>* pointers[5];
      std::vector<float> p(bins);
      for (int j = 0; j < T; ++j) {
        pointers[j] = inputs[j];
        for (int k = 0; k < bins; ++k) inputs[j][k] = {float(j+k+2), float(j-k)};
      }
      for (int k = 0; k < bins; ++k) p[k] = float(k+1)*3;
      NoisePower noise; noise.mode = PrimaryMode::Table; noise.table = {p.data(), p.size()}; noise.multiplier = float(T);
      select_fft3d_temporal(opt)(pointers, T, T/2, bins, 0, nullptr, noise, .1f, output);
      for (int k = 0; k < bins; ++k) {
        std::complex<double> result{};
        for (int m = 0; m < T; ++m) {
          std::complex<double> f{};
          for (int j = 0; j < T; ++j) f += std::complex<double>(inputs[j][k]) * std::polar(1., -2*3.141592653589793*j*m/T);
          const double q = std::norm(f)+1e-15;
          f *= std::max((q-T*p[k])/q, double(.1f));
          result += f*std::polar(1., 2*3.141592653589793*(T/2)*m/T)/double(T);
        }
        check_near(output[k].real(), result.real(), 2e-5);
        check_near(output[k].imag(), result.imag(), 2e-5);
      }
    }
    // Sampled power: independent direct DFT, rectangular windows, linear pfactor.
    for (float factor : {.5f, 1.f, 2.f}) {
      FFT3DConfig sampled; sampled.bw = sampled.bh = 4; sampled.ow = sampled.oh = 0;
      sampled.pfactor = factor; sampled.px = sampled.py = 1; sampled.degrid = 0; sampled.bt = 5;
      Plan plan(16, 16, {32,true,false}, sampled);
      std::vector<float> src(256);
      for (int y=0;y<16;++y) for(int x=0;x<16;++x) src[y*16+x] = float((x*3+y*7)%11)/8;
      const auto oracle = direct_dft(src.data(),4,4,16);
      const span2d::Plane<const float> plane{src.data(),16,16,64};
      auto one = std::async(std::launch::async,[&] { plan.prepare_pattern(plane); });
      auto two = std::async(std::launch::async,[&] { plan.prepare_pattern(plane); });
      one.get(); two.get();
      auto model = plan.pattern_power(); CHECK(model && model->size()==12);
      for(int y=0;y<4;++y) for(int x=0;x<3;++x) {
        const double r = std::pow(double(std::min(y,4-y))/2,2)+std::pow(double(x)/2,2);
        check_near((*model)[y*3+x], factor*std::norm(oracle[y*3+x])*r/(r+.01), 2e-4);
      }
      std::fill(src.begin(),src.end(),0); plan.prepare_pattern(plane);
      CHECK(plan.pattern_power()==model);
      CHECK(plan.workspace_pool().active_count()==0);
    }
    // Failure discards private state and releases the workspace; retry publishes.
    {
      FFT3DConfig sampled; sampled.bw=sampled.bh=4; sampled.ow=sampled.oh=0;
      sampled.pfactor=1; sampled.px=sampled.py=1;
      Plan plan(16,16,{32,true,false},sampled);
      std::vector<float> src(256,0); src[0]=INFINITY;
      span2d::Plane<const float> view{src.data(),16,16,64};
      rejects([&]{plan.prepare_pattern(view);}); CHECK(!plan.pattern_ready());
      CHECK(plan.workspace_pool().active_count()==0);
      src[0]=0; plan.prepare_pattern(view); CHECK(plan.pattern_ready());
      sampled.pshow=true;
      Plan float_preview(16,16,{32,true,false},sampled);
      std::vector<float> preview_output(256); src[255]=NAN;
      float_preview.process(view,{preview_output.data(),16,16,64});
      src[0]=NAN; rejects([&]{float_preview.process(view,{preview_output.data(),16,16,64});});
      sampled.px=sampled.py=0; sampled.pshow=true;
      Plan preview(16,16,{8,false,false},sampled);
      std::vector<std::uint8_t> input(256,128),output(256);
      preview.process({input.data(),16,16,16},{output.data(),16,16,16});
      for(int y=0;y<16;++y) for(int x=0;x<16;++x)
        CHECK(output[y*16+x]==(x>=4 && x<8 && y>=4 && y<8 ? 128 : 0));
      sampled.bw=sampled.bh=8; sampled.pshow=false;
      rejects([&]{Plan p(13,13,{8,false,false},sampled);});
      sampled.bt=-1; Plan inactive_sample(13,13,{8,false,false},sampled);
    }
    FFT3DConfig c; c.bt = -1; c.sigma = 1e30f;
    Plan identity(64, 64, {8, false, false}, c);
    c.bt = 1; rejects([&] { Plan p(64, 64, {8, false, false}, c); });
    c.bt = -1; c.enhancement.ht = 1e30f;
    Plan inactive(64, 64, {8, false, false}, c);
    c.enhancement.dehalo = 1; rejects([&] { Plan p(64, 64, {8, false, false}, c); });
    std::cout << "FFT3D profile, enhancement and temporal power oracles passed\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
