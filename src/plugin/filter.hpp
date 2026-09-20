#pragma once
#include "plugin/descriptors.hpp"
#include <dualsynth/video_filter.hpp>
#include <memory>
#include <cstring>

namespace neo_fft::plugin {
template <Algorithm A>
struct Filter {
  static constexpr const char* name = A == Algorithm::FFT3D ? "FFT3D" : "DFTTest";
  static constexpr int input_count = ds::dynamic_video_inputs;
  static constexpr ds::OutputOrigin output_origin = ds::OutputOrigin::fresh(0);
  struct State {
    ds::VideoInputInfo source;
    std::array<std::shared_ptr<const Plan>, 3> plans{};
  };
  static ds::Result<ds::VideoInitStateResult<State>> init(ds::VideoInitContext& ctx) {
    require(ctx.params && ctx.inputs.size() == 1, "missing clip or parameters");
    const auto& info = ctx.inputs[0];
    const auto& f = info.format;
    auto supported = ds::is_supported_video_format(f);
    require(supported.has_value() && supported.value() && (f.plane_count == 1 || f.plane_count == 3),
            "unsupported planar sample format");
    ds::validate_frame_dimensions(f, info.width, info.height);
    require(info.num_frames > 0, "clip frame count must be positive");
    Params params{*ctx.params};
    const auto config = [&] {
      if constexpr (A == Algorithm::FFT3D)
        return fft3d_config(params);
      else
        return dft_config(params);
    }();
    auto planes = params.integers("planes");
    std::array<bool, 3> selected{};
    if (!params.present("planes") || (A == Algorithm::FFT3D && planes.empty()))
      selected.fill(true);
    for (auto p : planes) {
      require(p >= 0 && p < f.plane_count, "planes index outside actual format");
      selected[std::size_t(p)] = true;
    }
    State state{info, {}};
    for (int p = 0; p < f.plane_count; ++p)
      if (selected[p]) {
        try {
          const bool chroma = f.color_family == ds::ColorFamily::Yuv && p > 0;
          const int w = info.width >> (chroma ? f.subsampling_w : 0), h = info.height >> (chroma ? f.subsampling_h : 0);
          state.plans[p] = std::make_shared<Plan>(
              w, h,
              SampleFormat{ds::bits_per_sample(f.sample_format), f.sample_format == ds::SampleFormat::Float32, chroma},
              config);
        } catch (const std::exception& e) {
          throw std::invalid_argument("plane " + std::to_string(p) + ": " + e.what());
        }
      }
    return ds::Result<ds::VideoInitStateResult<State>>::success(
        {{info.width, info.height, info.num_frames, f, info.fps}, std::move(state)});
  }
  static ds::VideoRequestPattern request_pattern(int, const State&) { return ds::VideoRequestPattern::StrictSpatial; }
  static ds::Result<ds::VideoRequestResult> request(ds::VideoRequestContext& ctx) {
    ctx.request_frame(0, ctx.output_frame);
    return ds::Result<ds::VideoRequestResult>::success({});
  }
  template <class T>
  static void plane(const ds::PlaneView& src, const ds::MutablePlaneView& dst, const Plan* plan) {
    const auto se = plane_extent<T>(src.width, src.height, src.stride_bytes);
    const auto de = plane_extent<T>(dst.width, dst.height, dst.stride_bytes);
    auto s = checked_plane(static_cast<const T*>(static_cast<const void*>(src.data)), src.width, src.height,
                           src.stride_bytes, se);
    auto d = checked_plane(static_cast<T*>(static_cast<void*>(dst.data)), dst.width, dst.height, dst.stride_bytes, de);
    disjoint(s.data(), se, d.data(), de);
    if (plan)
      plan->process(s, d);
    else
      for (int y = 0; y < src.height; ++y)
        std::memcpy(d.row_ptr(y), s.row_ptr(y), std::size_t(src.width) * sizeof(T));
  }
  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& ctx) {
    const auto& state = ctx.state<State>();
    auto source = unwrap(ctx.frames.get(0, ctx.output_frame));
    require(source.frame.format == state.source.format && ctx.dst.format == state.source.format &&
                source.frame.plane_count == state.source.format.plane_count &&
                ctx.dst.plane_count == source.frame.plane_count,
            "frame format differs from plan");
    for (int p = 0; p < source.frame.plane_count; ++p) {
      const auto& s = source.frame.plane(p);
      const auto& d = ctx.dst.plane(p);
      const bool chroma = state.source.format.color_family == ds::ColorFamily::Yuv && p > 0;
      const int w = state.source.width >> (chroma ? state.source.format.subsampling_w : 0);
      const int h = state.source.height >> (chroma ? state.source.format.subsampling_h : 0);
      require(s.width == w && s.height == h && d.width == w && d.height == h,
              "frame plane dimensions differ from plan");
      switch (state.source.format.sample_format) {
        case ds::SampleFormat::UInt8:
          plane<std::uint8_t>(s, d, state.plans[p].get());
          break;
        case ds::SampleFormat::Float32:
          plane<float>(s, d, state.plans[p].get());
          break;
        default:
          plane<std::uint16_t>(s, d, state.plans[p].get());
          break;
      }
    }
    return ds::Result<ds::VideoProcessResult>::success({});
  }
};
template <Algorithm A>
struct Bridge {
  using Core = Filter<A>;
  static constexpr const char* vs_name = A == Algorithm::FFT3D ? "FFT3D" : "DFTTest";
  inline static const std::string signature_storage = signature(A);
  inline static const char* vs_signature = signature_storage.c_str();
  static constexpr const char* avs_name = "";
  static constexpr const char* avs_signature = "";
  static constexpr const char* missing_input_error = "neo-fft: clip is required";
  static constexpr const char* vs_format_error = "neo-fft: fixed planar GRAY/YUV/RGB required";
  static constexpr const char* avs_format_error = vs_format_error;
  static constexpr std::size_t parity_source_index = 0;
  static constexpr bool forward_audio = false;
  static ds::FilterDescriptor descriptor() { return plugin::descriptor(A); }
  static bool accepts_video_format(const ds::VideoFormat& f) {
    auto v = ds::is_supported_video_format(f);
    return v.has_value() && v.value() && (f.plane_count == 1 || f.plane_count == 3);
  }
};
} // namespace neo_fft::plugin
