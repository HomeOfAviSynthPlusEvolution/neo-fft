#include "plugin/filter.hpp"
#include "plugin/avs_strings.hpp"
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
  static ds::FilterDescriptor descriptor(bool host_signature = false) {
    auto d = plugin::descriptor(A, host_signature);
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
    // Compatibility-only trailing arguments are deliberately not accessed.
    // Match neo-mv: native arrays occupy one slot; a scalar is a one-element array.
    std::vector<AVSValue> values(d.params.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = args.IsArray() ? (i < std::size_t(args.ArraySize()) ? args[int(i)] : AVSValue())
                                 : (i == 0 ? args : AVSValue());
      if constexpr (A == Algorithm::DFTTest) {
        if (accepts_dft_text_array(d.params[i].name) && values[i].IsString()) {
          std::vector<AVSValue> elements;
          const auto append = [&](const auto& numbers) {
            require(numbers.size() <= INT32_MAX, d.params[i].name + ": too many string elements");
            elements.reserve(numbers.size());
            for (const auto value : numbers)
              elements.emplace_back(value);
          };
          if (d.params[i].type == ds::ParamType::Integer)
            append(parse_number_list<int>(d.params[i].name.c_str(), values[i].AsString()));
          else
            append(parse_number_list<double>(d.params[i].name.c_str(), values[i].AsString()));
          values[i] = AVSValue(elements.data(), static_cast<int>(elements.size()));
        }
      }
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
  const auto signature = unwrap(ds::make_avisynth_signature(Adapter<A>::descriptor(true)));
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
