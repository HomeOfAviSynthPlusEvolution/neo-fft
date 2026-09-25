#pragma once
#include "algorithms/plan.hpp"
#include <dualsynth/param.hpp>

namespace neo_fft::plugin {
inline ds::FilterDescriptor descriptor(Algorithm algorithm, bool host_signature = false) {
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
    add("planes", P::Integer, true);
    for (auto n : {"bw", "bh", "bt", "ow", "oh"})
      add(n);
    for (auto n : {"kratio", "sharpen", "scutoff", "svr", "smin", "smax"})
      add(n, P::Float);
    add("interlaced", P::Boolean);
    for (auto n : {"wintype", "pframe", "px", "py"})
      add(n);
    add("pshow", P::Boolean);
    for (auto n : {"pcutoff", "pfactor", "sigma2", "sigma3", "sigma4", "degrid", "dehalo", "hr", "ht"})
      add(n, P::Float);
    for (auto n : {"l", "t", "r", "b", "opt"})
      add(n);
    add("cache_frames");
    add("cache_mb");
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
    add("opt");
  }
  // Signature-only compatibility arguments. Keep them after all consumed
  // arguments so neither host bridge parses them into algorithm parameters.
  if (host_signature) {
    if (algorithm == Algorithm::FFT3D) {
      add("mt", P::Boolean);
      add("ncpu");
      add("measure", P::Boolean);
    } else {
      add("threads");
      add("fft_threads");
    }
    add("fft_backend", P::String);
  }
  return d;
}
inline std::string signature(Algorithm a) {
  std::string out;
  for (const auto& p : descriptor(a, true).params) {
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
// Legacy modes enter ParamValues only through the AviSynth descriptor. An
// explicit modern selection takes precedence, including an empty array.
enum class PlaneMode { Skip = 1, Copy = 2, Process = 3 };
inline std::array<PlaneMode, 4> select_planes(Params params, Algorithm algorithm, int plane_count) {
  std::array<PlaneMode, 4> selected;
  selected.fill(PlaneMode::Copy);
  if (params.present("planes")) {
    const auto planes = params.integers("planes");
    if (algorithm == Algorithm::FFT3D && planes.empty())
      for (int p = 0; p < std::min(plane_count, 3); ++p) selected[p] = PlaneMode::Process;
    for (auto p : planes) {
      require(p >= 0 && p < plane_count, "planes index outside actual format");
      selected[std::size_t(p)] = PlaneMode::Process;
    }
  } else {
    constexpr const char* names[] = {"y", "u", "v", "a"};
    for (int p = 0; p < 4; ++p) {
      const int mode = params.integer(names[p], p == 3 ? 2 : 3);
      require(mode >= 1 && mode <= 3, std::string(names[p]) + ": mode must be 1, 2 or 3");
      if (p < plane_count) selected[p] = static_cast<PlaneMode>(mode);
    }
  }
  return selected;
}
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

  c.kratio=p.number("kratio",2.0f);
  c.enhancement.sharpen = p.number("sharpen", 0.0f);
  c.enhancement.scutoff = p.number("scutoff", 0.3f);
  c.enhancement.svr = p.number("svr", 1.0f);
  c.enhancement.smin = p.number("smin", 4.0f);
  c.enhancement.smax = p.number("smax", 20.0f);
  c.interlaced = p.boolean("interlaced", false);
  c.left = p.integer("l", 0); c.top = p.integer("t", 0);
  c.right = p.integer("r", 0); c.bottom = p.integer("b", 0);
  c.pframe = p.integer("pframe", 0);
  c.px = p.integer("px", 0); c.py = p.integer("py", 0);
  c.pshow = p.boolean("pshow", false);
  c.pcutoff = p.number("pcutoff", .1f);
  c.pfactor = p.number("pfactor", 0);
  c.sigma2 = p.number("sigma2", c.sigma);
  c.sigma3 = p.number("sigma3", c.sigma);
  c.sigma4 = p.number("sigma4", c.sigma);
  c.enhancement.dehalo = p.number("dehalo", 0.0f);
  c.enhancement.hr = p.number("hr", 2.0f);
  c.enhancement.ht = p.number("ht", 50.0f);
  c.cache_frames=p.integer("cache_frames",-1);
  c.cache_mb=p.integer("cache_mb",runtime::SpectraCache::default_mb);
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
  c.temporal_mode = p.integer("tmode",0);
  const int temporal_overlap = p.integer("tosize",0);
  c.temporal_overlap = c.temporal_mode==0 ? 0 : temporal_overlap;
  const auto locations = p.integers("nlocation");
  require(locations.size()%4 == 0 && locations.size()/4 <= 500, "DFTTest nlocation requires up to 500 quadruples");
  for (std::size_t i=0;i<locations.size();i+=4)
    c.locations.push_back({int(locations[i]),int(locations[i+1]),int(locations[i+2]),int(locations[i+3])});
  c.alpha = p.number("alpha", c.ftype == 0 ? 5.0f : 7.0f);
  c.curves.shared = p.numbers("slocation");
  c.curves.x = p.numbers("ssx");
  c.curves.y = p.numbers("ssy");
  c.curves.time = p.numbers("sst");
  c.curves.system = p.integer("ssystem", 0);
  c.dither=p.integer("dither",0);
  c.dither_seed=p.integer("dither_seed",0);
  validate(c);
  return c;
}
} // namespace neo_fft::plugin
