#include "algorithms/profile.hpp"
#include <algorithm>
#include "kernels/table.hpp"

namespace neo_fft {
namespace {
using Curve = std::vector<std::pair<float, float>>;
Curve sorted(const std::vector<float>& values) {
  require(values.size() % 2 == 0, "DFTTest curve requires frequency/value pairs");
  Curve knots;
  for (std::size_t i = 0; i < values.size(); i += 2) {
    require(std::isfinite(values[i]) && values[i] >= 0 && values[i] <= 1, "DFTTest curve position outside 0..1");
    require(std::isfinite(values[i + 1]) && values[i + 1] >= 0, "DFTTest invalid curve value");
    knots.emplace_back(values[i], values[i + 1]);
  }
  std::sort(knots.begin(), knots.end());
  if (!knots.empty()) {
    require(knots.front().first == 0 && knots.back().first == 1, "DFTTest curve requires endpoints 0 and 1");
    for (std::size_t i = 1; i < knots.size(); ++i)
      require(knots[i - 1].first < knots[i].first, "DFTTest duplicate curve position after binary32 conversion");
  }
  return knots;
}
Curve prepare(const std::vector<float>& values, float sigma, float exponent) {
  auto knots = sorted(values);
  if (knots.empty())
    knots = {{0, sigma}, {1, sigma}};
  for (auto& knot : knots)
    knot.second = finite(std::pow(knot.second, exponent));
  return knots;
}
float interpolate(const Curve& curve, float f) {
  const auto hi =
      std::lower_bound(curve.begin(), curve.end(), f, [](const auto& knot, float value) { return knot.first < value; });
  require(hi != curve.end(), "DFTTest frequency outside curve");
  if (hi->first == f)
    return hi->second;
  require(hi != curve.begin(), "DFTTest frequency outside curve");
  const auto lo = hi - 1;
  const float t = finite((f - lo->first) / (hi->first - lo->first));
  return finite(finite(lo->second * (1 - t)) + finite(hi->second * t));
}
float frequency(int i, int length) {
  return length == 1 ? 0 : float(std::min(i, length - i)) / float(length / 2);
}
} // namespace
void validate(const DFTCurves& c) {
  require(c.system == 0 || c.system == 1, "DFTTest ssystem outside 0..1");
  for (const auto* values : {&c.shared, &c.x, &c.y, &c.time})
    sorted(*values);
}
std::vector<float> dft_profile(const DFTCurves& c, int T, int S, float sigma, float divisor, int opt) {
  validate(c);
  require(T > 0 && S > 0 && std::isfinite(sigma) && sigma >= 0 && std::isfinite(divisor) && divisor > 0,
          "DFTTest invalid profile shape or scale");
  auto result = buffer<float>(mul_size(mul_size(std::size_t(T), std::size_t(S)), std::size_t(S / 2 + 1)));
  if (c.empty()) {
    std::fill(result.begin(), result.end(), finite(sigma / divisor));
    return result;
  }
  const int dimensions = (T > 1 ? 1 : 0) + (S > 1 ? 2 : 0);
  const bool shared = !c.shared.empty();
  const float exponent = dimensions == 0 || (shared && c.system == 1) ? 1.0f : 1.0f / float(dimensions);
  const auto time = prepare(shared ? c.shared : c.time, sigma, exponent);
  // The all-singleton extension consumes only the temporal/shared DC value.
  if (dimensions == 0) {
    result[0] = finite(interpolate(time, 0) / divisor);
    return result;
  }
  const auto y = prepare(shared ? c.shared : c.y, sigma, exponent);
  const auto x = prepare(shared ? c.shared : c.x, sigma, exponent);
  const auto kernels=select_table(opt);
  const int K=S/2+1;
  auto ft=buffer<float>(T),fy=buffer<float>(S),fx=buffer<float>(K);
  for(int z=0;z<T;++z)ft[z]=frequency(z,T);
  for(int j=0;j<S;++j)fy[j]=frequency(j,S);
  for(int i=0;i<K;++i)fx[i]=frequency(i,S);
  auto evaluate=[&](const Curve& knots,const std::vector<float>& freq,std::vector<float>& out,float scale) {
    auto positions=buffer<float>(knots.size()),values=buffer<float>(knots.size());
    for(std::size_t i=0;i<knots.size();++i) {positions[i]=knots[i].first;values[i]=knots[i].second;}
    kernels.curve(out.data(),freq.data(),freq.size(),positions.data(),values.data(),knots.size(),scale);
  };
  if(c.system==1) {
    auto radius=buffer<float>(K);
    auto positions=buffer<float>(time.size()),values=buffer<float>(time.size());
    for(std::size_t i=0;i<time.size();++i) {positions[i]=time[i].first;values[i]=time[i].second;}
    for(int z=0;z<T;++z) for(int j=0;j<S;++j) {
      kernels.radius(radius.data(),fx.data(),fx.size(),ft[z]*ft[z]+fy[j]*fy[j],float(dimensions));
      kernels.curve(result.data()+(std::size_t(z)*S+j)*K,radius.data(),radius.size(),
                    positions.data(),values.data(),time.size(),divisor);
    }
  } else {
    auto vt=buffer<float>(T),vy=buffer<float>(S),vx=buffer<float>(K);
    if(T==1)vt[0]=1;else evaluate(time,ft,vt,1);
    if(S==1) {vy[0]=1;vx[0]=1;} else {evaluate(y,fy,vy,1);evaluate(x,fx,vx,1);}
    for(int z=0;z<T;++z) for(int j=0;j<S;++j)
      kernels.product(result.data()+(std::size_t(z)*S+j)*K,vx.data(),vx.size(),finite(vt[z]*vy[j]),divisor);
  }
  return result;
}
} // namespace neo_fft
