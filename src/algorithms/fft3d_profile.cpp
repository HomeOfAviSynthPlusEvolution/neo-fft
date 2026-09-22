#include "algorithms/fft3d_profile.hpp"
#include <algorithm>
#include "kernels/table.hpp"

namespace neo_fft {
void validate(const EnhancementConfig& c) {
  for (float v : {c.sharpen, c.svr, c.smin, c.smax, c.dehalo, c.ht})
    require(std::isfinite(v) && v >= 0, "FFT3D invalid enhancement parameter");
  require(c.smin <= c.smax, "FFT3D smin exceeds smax");
  for (float v : {c.scutoff, c.hr})
    require(std::isfinite(v) && v > 0, "FFT3D cutoff/radius must be positive");
}
std::vector<float> fft3d_profile(int W, int H, const std::array<float, 4>& s, int opt) {
  require(W >= 2 && H >= 2, "FFT3D profile shape invalid");
  for (float v : s) require(std::isfinite(v) && v >= 0, "FFT3D invalid scaled sigma");
  const int K = W / 2 + 1;
  const float norm = 1 / (float(W) * float(H));
  const auto kernels=select_table(opt);
  auto p = buffer<float>(mul_size(std::size_t(H), std::size_t(K)));
  auto fx=buffer<float>(K),radius=buffer<float>(K);
  for(int x=0;x<K;++x) fx[x]=float(x)/float(K);
  for(int y=0;y<H;++y) {
    const float fy=float(H-2*std::abs(y-H/2))/float(H);
    kernels.radius(radius.data(),fx.data(),fx.size(),fy*fy,2);
    kernels.analytic(p.data()+std::size_t(y)*K,radius.data(),radius.size(),s.data(),norm);
  }
  return p;
}
Enhancement EnhancementTables::view(const EnhancementConfig& config) const {
  return {config.sharpen, config.dehalo, a, b, c, {sharpen.data(), sharpen.size()}, {halo.data(), halo.size()}};
}
EnhancementTables enhancement_tables(int W, int H, float factor, const EnhancementConfig& c, int opt) {
  validate(c);
  require(W >= 2 && H >= 2 && std::isfinite(factor) && factor > 0, "FFT3D invalid enhancement geometry");
  EnhancementTables result;
  if (c.sharpen == 0 && c.dehalo == 0) return result;
  const float norm = 1 / (float(W) * float(H));
  float cutoff = 1;
  if (c.sharpen != 0) {
    const float lo = finite(c.smin * factor), hi = finite(c.smax * factor);
    result.a = finite(finite(lo * lo) / norm);
    result.b = finite(finite(hi * hi) / norm);
    cutoff = finite(finite(2 * c.scutoff) * c.scutoff);
    require(cutoff > 0, "FFT3D sharpen cutoff underflow");
  }
  if (c.dehalo != 0) result.c = finite(finite(c.ht * c.ht) / norm);
  const auto count = mul_size(std::size_t(H), std::size_t(W / 2 + 1));
  if (c.sharpen != 0) result.sharpen = buffer<float>(count);
  if (c.dehalo != 0) result.halo = buffer<float>(count);
  const auto kernels=select_table(opt);
  const int K=W/2+1;
  auto x2=buffer<float>(K);
  for(int x=0;x<K;++x) x2[x]=float(x)*float(x)/(float(W/2)*float(W/2));
  for(int y=0;y<H;++y) {
    const float dy=float(y<H/2 ? y : H-y);
    const float y2=finite(finite(finite(dy*dy)*c.svr)*c.svr)/(float(H/2)*float(H/2));
    const auto offset=std::size_t(y)*K;
    kernels.enhancement(c.sharpen!=0 ? result.sharpen.data()+offset : nullptr,
                        c.dehalo!=0 ? result.halo.data()+offset : nullptr,
                        x2.data(),x2.size(),y2,cutoff,c.hr);
  }
  const float maximum=kernels.maximum(result.halo.data(),result.halo.size());
  if (c.dehalo != 0) {
    require(std::isfinite(maximum) && maximum > 0, "FFT3D degenerate dehalo window");
    kernels.divide(result.halo.data(),result.halo.size(),maximum);
  }
  return result;
}
float enhancement_gain(float q, std::size_t k, const Enhancement& e) {
  float s = 1, d = 1;
  if (e.sharpen != 0 && e.b != 0) {
    const float num = finite(q * e.b), den = finite(finite(q + e.a) * finite(q + e.b));
    const float strength = finite(e.sharpen * e.sharpen_window[k]);
    s = finite(1 + finite(strength * finite(std::sqrt(finite(num / den)))));
  }
  if (e.dehalo != 0) {
    const float qc = finite(q + e.c);
    d = finite(qc / finite(qc + finite(finite(e.dehalo * e.halo_window[k]) * q)));
  }
  return finite(s * d);
}
} // namespace neo_fft
