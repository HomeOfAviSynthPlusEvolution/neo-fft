#pragma once
#include "algorithms/plan.hpp"
#include <dualsynth/param.hpp>

namespace neo_fft::plugin {
inline ds::FilterDescriptor descriptor(Algorithm algorithm) {
  using P = ds::ParamType;
  ds::FilterDescriptor d;
  d.name = algorithm == Algorithm::FFT3D ? "FFT3D" : "DFTTest";
  auto add = [&](const char* name, P type = P::Integer, bool array = false) {
    d.params.push_back({name, type, {}, std::string(name) == "clip", array, true, false});
  };
  add("clip", P::Clip);
  if (algorithm == Algorithm::FFT3D) {
    add("sigma", P::Float);
    add("beta", P::Float);
    add("fft_backend", P::String);
    add("planes", P::Integer, true);
    for (auto n : {"bw", "bh", "bt", "ow", "oh"})
      add(n);
    for (auto n : {"kratio", "sharpen", "scutoff", "svr", "smin", "smax"})
      add(n, P::Float);
    add("measure", P::Boolean);
    add("interlaced", P::Boolean);
    for (auto n : {"wintype", "pframe", "px", "py"})
      add(n);
    add("pshow", P::Boolean);
    for (auto n : {"pcutoff", "pfactor", "sigma2", "sigma3", "sigma4", "degrid", "dehalo", "hr", "ht"})
      add(n, P::Float);
    for (auto n : {"l", "t", "r", "b", "opt", "ncpu"})
      add(n);
    add("mt", P::Boolean);
  } else {
    add("ftype");
    for (auto n : {"sigma", "sigma2", "pmin", "pmax"})
      add(n, P::Float);
    for (auto n : {"sbsize", "smode", "sosize", "tbsize", "tmode", "tosize", "swin", "twin"})
      add(n);
    add("sbeta", P::Float);
    add("tbeta", P::Float);
    add("zmean", P::Boolean);
    add("f0beta", P::Float);
    add("nlocation", P::Integer, true);
    add("alpha", P::Float);
    for (auto n : {"slocation", "ssx", "ssy", "sst"})
      add(n, P::Float, true);
    for (auto n : {"ssystem", "dither", "dither_seed"})
      add(n);
    add("planes", P::Integer, true);
    for (auto n : {"opt", "threads", "fft_threads"})
      add(n);
    add("fft_backend", P::String);
  }
  return d;
}
inline std::string signature(Algorithm a) {
  std::string out;
  for (const auto& p : descriptor(a).params) {
    out += p.name + ":";
    switch (p.type) {
      case ds::ParamType::Clip:
        out += "vnode";
        break;
      case ds::ParamType::Float:
        out += "float";
        break;
      case ds::ParamType::Integer:
      case ds::ParamType::Boolean:
        out += "int";
        break;
      case ds::ParamType::String:
        out += "data";
        break;
    }
    if (p.is_array)
      out += "[]";
    if (!p.required)
      out += ":opt";
    if (p.is_array)
      out += ":empty";
    out += ";";
  }
  return out;
}
template <class T>
T unwrap(ds::Result<T> result) {
  if (!result.has_value())
    throw std::invalid_argument(result.error().message);
  return std::move(result.value());
}
struct Params {
  const ds::ParamValues& values;
  bool present(const char* name) const {
    for (const auto& v : values.entries)
      if (v.name == name)
        return true;
    return false;
  }
  int integer(const char* n, int d) const {
    auto v = unwrap(values.get_int64(n, d));
    require(v >= INT32_MIN && v <= INT32_MAX, std::string(n) + ": integer outside int32");
    return static_cast<int>(v);
  }
  float number(const char* n, float d) const {
    const double v = unwrap(values.get_double(n, d));
    require(std::isfinite(v) && std::abs(v) <= std::numeric_limits<float>::max(),
            std::string(n) + ": non-finite or unrepresentable float");
    return float(v);
  }
  bool boolean(const char* n, bool d) const { return unwrap(values.get_bool(n, d)); }
  void same(const char* n, int d) const { require(integer(n, d) == d, std::string("unsupported phase-1 ") + n); }
  void same(const char* n, float d) const { require(number(n, d) == d, std::string("unsupported phase-1 ") + n); }
  void off(const char* n) const { require(!boolean(n, false), std::string("unsupported phase-1 ") + n); }
  void backend() const {
    require(unwrap(values.get_string("fft_backend", "pocketfft")) == "pocketfft",
            "unsupported fft_backend: only pocketfft is built");
  }
  std::vector<std::int64_t> integers(const char* n) const {
    auto v = unwrap(values.get_int_array(n, {}));
    for (auto x : v)
      require(x >= INT32_MIN && x <= INT32_MAX, std::string(n) + ": array element outside int32");
    return v;
  }
  std::vector<float> numbers(const char* n) const {
    auto v = unwrap(values.get_double_array(n, {}));
    for (auto x : v)
      require(std::isfinite(x) && std::abs(x) <= std::numeric_limits<float>::max(),
              std::string(n) + ": invalid array float");
    return std::vector<float>(v.begin(), v.end());
  }
};
inline FFT3DConfig fft3d_config(Params p) {
  FFT3DConfig c;
  c.bw = p.integer("bw", 32);
  c.bh = p.integer("bh", 32);
  c.ow = p.integer("ow", -1);
  c.oh = p.integer("oh", -1);
  c.wintype = p.integer("wintype", 0);
  c.opt = p.integer("opt", 0);
  c.sigma = p.number("sigma", 2);
  c.beta = p.number("beta", 1);
  c.degrid = p.number("degrid", 1);
  c.bt = p.integer("bt", 3);
  require(c.bt >= 1 && c.bt <= 5, "FFT3D bt outside 1..5");
  p.same("kratio", 2.0f);
  p.same("sharpen", 0.0f);
  p.same("scutoff", 0.3f);
  p.same("svr", 1.0f);
  p.same("smin", 4.0f);
  p.same("smax", 20.0f);
  p.boolean("measure", true);
  p.off("interlaced");
  for (auto n : {"pframe", "px", "py", "l", "t", "r", "b"})
    p.same(n, 0);
  p.off("pshow");
  p.same("pcutoff", 0.1f);
  p.same("pfactor", 0.0f);
  for (auto n : {"sigma2", "sigma3", "sigma4"})
    p.same(n, c.sigma);
  p.same("dehalo", 0.0f);
  p.same("hr", 2.0f);
  p.same("ht", 50.0f);
  require(p.integer("ncpu", 2) > 0, "ncpu must be positive");
  p.off("mt");
  p.backend();
  validate(c);
  return c;
}
inline DFTConfig dft_config(Params p) {
  DFTConfig c;
  c.ftype = p.integer("ftype", 0);
  c.sigma = p.number("sigma", 8);
  c.sigma2 = p.number("sigma2", 8);
  c.pmin = p.number("pmin", 0);
  c.pmax = p.number("pmax", 500);
  c.block = p.integer("sbsize", 16);
  c.mode = p.integer("smode", 1);
  c.overlap = p.integer("sosize", 12);
  c.swin = p.integer("swin", 0);
  c.twin = p.integer("twin", 7);
  c.opt = p.integer("opt", 0);
  c.sbeta = p.number("sbeta", 2.5f);
  c.tbeta = p.number("tbeta", 2.5f);
  c.f0beta = p.number("f0beta", 1);
  c.zmean = p.boolean("zmean", true);
  c.tbsize = p.integer("tbsize", 3);
  require(c.tbsize >= 1 && c.tbsize <= 15 && c.tbsize % 2 == 1, "DFTTest tbsize must be odd integer in 1..15");
  p.same("tmode", 0);
  p.integer("tosize", 0);
  require(p.integers("nlocation").empty(), "unsupported phase-1 nlocation");
  p.same("alpha", c.ftype == 0 ? 5.0f : 7.0f);
  c.curves.shared = p.numbers("slocation");
  c.curves.x = p.numbers("ssx");
  c.curves.y = p.numbers("ssy");
  c.curves.time = p.numbers("sst");
  c.curves.system = p.integer("ssystem", 0);
  p.same("dither", 0);
  require(p.integer("dither_seed", 0) >= 0, "dither_seed must be nonnegative");
  require(p.integer("threads", 0) <= 1, "unsupported phase-1 threads > 1");
  require(p.integer("fft_threads", 0) <= 1, "unsupported phase-1 fft_threads > 1");
  p.backend();
  validate(c);
  return c;
}
} // namespace neo_fft::plugin
