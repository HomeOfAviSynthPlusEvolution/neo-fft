#include "plugin/filter.hpp"
#include <avisynth.h>
#include <dualsynth/avisynth/video_bridge.hpp>

const AVS_Linkage* AVS_linkage = nullptr;

namespace neo_fft::plugin::avs {
namespace av = ds::avisynth;

template <Algorithm A>
struct Adapter : Bridge<A> {
  static constexpr const char* avs_name = A == Algorithm::FFT3D ? "neo_fft_FFT3D" : "neo_fft_DFTTest";
  static constexpr bool forward_audio = true;
  static constexpr av::MtMode avs_mt_mode = av::MtMode::NiceFilter;
  static ds::FilterDescriptor descriptor() {
    auto d = Bridge<A>::descriptor();
    for (auto& p : d.params) {
      p.avs_enabled = true;
      if (p.is_array)
        p.avs_array_binding = ds::AvisynthArrayBinding::Native;
    }
    return d;
  }
};

template <class Function>
AVSValue guarded(IScriptEnvironment* env, Function&& function) {
  try {
    env->CheckVersion(11);
    return function();
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo-fft: %s", e.what());
  } catch (...) {
    env->ThrowError("neo-fft: AviSynth call failed");
  }
  return {};
}

template <Algorithm A>
AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env) {
  return guarded(env, [&] {
    const auto d = Adapter<A>::descriptor();
    // Match neo-mv: native arrays occupy one slot; a scalar is a one-element array.
    std::vector<AVSValue> values(d.params.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = args.IsArray() ? (i < std::size_t(args.ArraySize()) ? args[int(i)] : AVSValue())
                                 : (i == 0 ? args : AVSValue());
      if (d.params[i].is_array && values[i].Defined() && !values[i].IsArray()) {
        const AVSValue scalar = values[i];
        values[i] = AVSValue(&scalar, 1);
      }
    }
    return av::create_video_filter_bridge<Adapter<A>>(
        AVSValue(values.data(), static_cast<int>(values.size())), env);
  });
}

AVSValue __cdecl kernel_info(AVSValue, void*, IScriptEnvironment* env) {
  return guarded(env, [&] {
    AVSValue fields[] = {AVSValue(env->SaveString(fft_backend_name())),
                         AVSValue(env->SaveString(spectral_target(0))),
                         AVSValue(env->SaveString(fft_profile_name())), AVSValue(fft_lanes())};
    return AVSValue(fields, 4);
  });
}

template <Algorithm A>
void add(IScriptEnvironment* env) {
  const auto signature = unwrap(ds::make_avisynth_signature(Adapter<A>::descriptor()));
  env->AddFunction(Adapter<A>::avs_name, env->SaveString(signature.c_str()), create<A>, nullptr);
}
} // namespace neo_fft::plugin::avs

#if defined(_WIN32)
#define NEO_FFT_AVS_EXPORT extern "C" __declspec(dllexport)
#else
#define NEO_FFT_AVS_EXPORT extern "C" __attribute__((visibility("default")))
#endif

NEO_FFT_AVS_EXPORT const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* linkage) {
  AVS_linkage = linkage;
  neo_fft::plugin::avs::guarded(env, [&] {
    neo_fft::plugin::avs::add<neo_fft::Algorithm::FFT3D>(env);
    neo_fft::plugin::avs::add<neo_fft::Algorithm::DFTTest>(env);
    env->AddFunction("neo_fft_KernelInfo", "", neo_fft::plugin::avs::kernel_info, nullptr);
    return AVSValue();
  });
  return "neo-fft AviSynth interface";
}
