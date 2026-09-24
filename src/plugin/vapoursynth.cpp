#include "plugin/filter.hpp"
#include <dualsynth/vapoursynth/video_bridge.hpp>

namespace neo_fft::plugin {
namespace {
// Check consumed parameters, including inactive arrays. Signature-only
// compatibility arguments are absent here and from the DS2 bridge descriptor.
void validate_supplied(const VSMap* in, const VSAPI* api, Algorithm a) {
  for (const auto& p : descriptor(a).params) {
    const int count = api->mapNumElements(in, p.name.c_str());
    for (int i = 0; i < count; ++i) {
      int error = 0;
      if (p.type == ds::ParamType::Integer) {
        const auto v = api->mapGetInt(in, p.name.c_str(), i, &error);
        require(!error && v >= INT32_MIN && v <= INT32_MAX, p.name + ": integer outside int32");
      } else if (p.type == ds::ParamType::Float) {
        const double v = api->mapGetFloat(in, p.name.c_str(), i, &error);
        require(!error && std::isfinite(v) && std::abs(v) <= std::numeric_limits<float>::max(),
                p.name + ": invalid float");
      } else if (p.type == ds::ParamType::String) {
        const auto* data = api->mapGetData(in, p.name.c_str(), i, &error);
        const int size = api->mapGetDataSize(in, p.name.c_str(), i, &error);
        require(!error && data && size >= 0 && !std::memchr(data, 0, std::size_t(size)), p.name + ": invalid string");
      }
    }
  }
}
template <Algorithm A>
void VS_CC create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* api) {
  try {
    validate_supplied(in, api, A);
    if constexpr(A==Algorithm::FFT3D) {
      // VS has no AviSynth cache-hint notification. Resolve auto capacity from
      // the core's configured worker count at creation, without mutating input.
      const int supplied=api->mapNumElements(in,"cache_frames");
      require(supplied==-1 || supplied==1,"cache_frames: expected one integer");
      int error=0;
      const auto frames=api->mapGetInt(in,"cache_frames",0,&error);
      if(error || frames==-1) {
        VSCoreInfo info{};api->getCoreInfo(core,&info);
        const auto bt_value=api->mapGetInt(in,"bt",0,&error);
        const int bt=error ? 3 : int(bt_value);
        const auto count=std::min<std::int64_t>(INT32_MAX,std::int64_t(std::max(1,bt))+std::max(1,info.numThreads)-1);
        const auto free_map=[&](VSMap* map){api->freeMap(map);};
        std::unique_ptr<VSMap,decltype(free_map)> params(api->createMap(),free_map);
        require(bool(params),"cache auto parameters allocation failed");
        api->copyMap(in,params.get());
        require(!api->mapSetInt(params.get(),"cache_frames",count,maReplace),"cache auto parameter failed");
        ds::vapoursynth::create_video_filter_bridge<Bridge<A>>(params.get(),out,core,api);
      } else ds::vapoursynth::create_video_filter_bridge<Bridge<A>>(in,out,core,api);
    } else ds::vapoursynth::create_video_filter_bridge<Bridge<A>>(in, out, core, api);
    if (const auto* e = api->mapGetError(out)) {
      const std::string owned = std::string(Bridge<A>::vs_name) + ": " + e;
      api->mapSetError(out, owned.c_str());
    }
  } catch (const std::exception& e) {
    const std::string owned = std::string(Bridge<A>::vs_name) + ": " + e.what();
    api->mapSetError(out, owned.c_str());
  } catch (...) {
    api->mapSetError(out, "neo-fft: unhandled creation failure");
  }
}
void VS_CC info(const VSMap*, VSMap* out, void*, VSCore*, const VSAPI* api) {
  api->mapSetData(out, "fft_backend", fft_backend_name(), -1, dtUtf8, maReplace);
  api->mapSetData(out, "target", spectral_target(0), -1, dtUtf8, maReplace);
  api->mapSetData(out, "fft", fft_profile_name(), -1, dtUtf8, maReplace);
  api->mapSetInt(out, "fft_lanes", fft_lanes(), maReplace);
}
} // namespace
} // namespace neo_fft::plugin
VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* api) {
  using namespace neo_fft;
  using namespace neo_fft::plugin;
  api->configPlugin("org.neofilters.neo_fft", "neo_fft", "neo-fft spatial filters", VS_MAKE_VERSION(0, 1),
                    VAPOURSYNTH_API_VERSION, 0, plugin);
  const auto f3d = signature(Algorithm::FFT3D), dft = signature(Algorithm::DFTTest);
  api->registerFunction("FFT3D", f3d.c_str(), "clip:vnode;", create<Algorithm::FFT3D>, nullptr, plugin);
  api->registerFunction("DFTTest", dft.c_str(), "clip:vnode;", create<Algorithm::DFTTest>, nullptr, plugin);
  api->registerFunction("KernelInfo", "",
                        "fft_backend:data;target:data;fft:data:opt;fft_lanes:int:opt;", info, nullptr,
                        plugin);
}
