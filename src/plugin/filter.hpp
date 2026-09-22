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
    int temporal_size = 1;
    int pattern_frame = 0;
    bool sampled = false;
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
    require(info.num_frames <= std::numeric_limits<int>::max() - 16, "clip frame count exceeds supported limit");
    Params params{*ctx.params};
    const auto config = [&] {
      if constexpr (A == Algorithm::FFT3D)
        return fft3d_config(params);
      else
        return dft_config(params);
    }();
    if constexpr (A == Algorithm::DFTTest) {
      require(config.tbsize <= info.num_frames, "tbsize must be less than or equal to the number of frames");
    }
    auto planes = params.integers("planes");
    std::array<bool, 3> selected{};
    if (!params.present("planes") || (A == Algorithm::FFT3D && planes.empty()))
      selected.fill(true);
    for (auto p : planes) {
      require(p >= 0 && p < f.plane_count, "planes index outside actual format");
      selected[std::size_t(p)] = true;
    }
    const int t_size = [&] {
      if constexpr (A == Algorithm::FFT3D)
        return std::max(1, config.bt);
      else
        return config.tbsize;
    }();
    State state{info, {}, t_size};
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
    if constexpr (A == Algorithm::FFT3D) {
      state.pattern_frame = std::clamp(config.pframe, 0, info.num_frames-1);
      for (const auto& plan : state.plans) if (plan) {
        if (plan->preview()) state.temporal_size = 1;
        state.sampled = state.sampled || plan->needs_pattern_frame();
      }
    }
    return ds::Result<ds::VideoInitStateResult<State>>::success(
        {{info.width, info.height, info.num_frames, f, info.fps}, std::move(state)});
  }
  static ds::VideoRequestPattern request_pattern(int, const State& state) {
    return state.temporal_size > 1 || state.sampled ? ds::VideoRequestPattern::General : ds::VideoRequestPattern::StrictSpatial;
  }
  static ds::Result<ds::VideoRequestResult> request(ds::VideoRequestContext& ctx) {
    const auto& state = ctx.state<State>();
    const int n = ctx.output_frame;
    const int N = state.source.num_frames;

    const bool has_active_plans = std::any_of(state.plans.begin(), state.plans.end(), [](const auto& p) {
      return p != nullptr;
    });
    if (!has_active_plans) {
      ctx.request_frame(0, n);
      return ds::Result<ds::VideoRequestResult>::success({});
    }

    if constexpr (A == Algorithm::FFT3D) {
      if (state.sampled) ctx.request_frame(0, state.pattern_frame);
      const int bt = state.temporal_size;
      const int left = bt / 2;
      const int right = (bt - 1) / 2;
      const bool fallback = (bt <= 1) || (n < left || N - 1 - n < right);
      if (fallback) {
        ctx.request_frame(0, n);
      } else {
        for (int f = n - left; f <= n + right; ++f) {
          ctx.request_frame(0, f);
        }
      }
    } else {
      const int T = state.temporal_size;
      const int c = T / 2;
      for (int j = 0; j < T; ++j) {
        const int real = (j >= c) ? ((N - 1 - n < j - c) ? (N - 1) : (n + (j - c)))
                                   : ((n < c - j) ? 0 : (n - (c - j)));
        ctx.request_frame(0, real);
      }
    }
    return ds::Result<ds::VideoRequestResult>::success({});
  }
  template <class T>
  static void plane(span2d::Span<const ds::PlaneView> src_views, const ds::MutablePlaneView& dst, const Plan* plan) {
    const auto de = plane_extent<T>(dst.width, dst.height, dst.stride_bytes);
    auto d = checked_plane(static_cast<T*>(static_cast<void*>(dst.data)), dst.width, dst.height, dst.stride_bytes, de);

    if (plan) {
      std::vector<span2d::Plane<const T>> checked_srcs;
      checked_srcs.reserve(src_views.size());
      for (const auto& src : src_views) {
        const auto se = plane_extent<T>(src.width, src.height, src.stride_bytes);
        checked_srcs.push_back(checked_plane(static_cast<const T*>(static_cast<const void*>(src.data)),
                                             src.width, src.height, src.stride_bytes, se));
        disjoint(checked_srcs.back().data(), se, d.data(), de);
      }
      plan->process(span2d::Span<const span2d::Plane<const T>>(checked_srcs.data(), checked_srcs.size()), d);
    } else {
      const int c = int(src_views.size()) / 2;
      const auto& src = src_views[c];
      const auto se = plane_extent<T>(src.width, src.height, src.stride_bytes);
      const auto s = checked_plane(static_cast<const T*>(static_cast<const void*>(src.data)),
                                   src.width, src.height, src.stride_bytes, se);
      disjoint(s.data(), se, d.data(), de);
      for (int y = 0; y < src.height; ++y)
        std::memcpy(d.row_ptr(y), s.row_ptr(y), std::size_t(src.width) * sizeof(T));
    }
  }
  template<class T> static void prepare_pattern(const ds::PlaneView& src, const Plan& plan) {
    const auto s = checked_plane(static_cast<const T*>(static_cast<const void*>(src.data)), src.width, src.height,
        src.stride_bytes, plane_extent<T>(src.width, src.height, src.stride_bytes));
    plan.prepare_pattern(s);
  }
  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& ctx) {
    const auto& state = ctx.state<State>();
    const int n = ctx.output_frame;
    const int N = state.source.num_frames;

    const bool has_active_plans = std::any_of(state.plans.begin(), state.plans.end(), [](const auto& p) {
      return p != nullptr;
    });

    int T = 1;
    std::vector<ds::RequestedVideoFrame> frames_holder;

    if (!has_active_plans) {
      frames_holder.reserve(1);
      frames_holder.push_back(unwrap(ctx.frames.get(0, n)));
    } else if constexpr (A == Algorithm::FFT3D) {
      const int bt = state.temporal_size;
      const int left = bt / 2;
      const int right = (bt - 1) / 2;
      const bool fallback = (bt <= 1) || (n < left || N - 1 - n < right);
      T = fallback ? 1 : bt;
      const int c = T / 2;
      frames_holder.reserve(T);
      for (int j = 0; j < T; ++j) {
        const int f = fallback ? n : (n - c + j);
        frames_holder.push_back(unwrap(ctx.frames.get(0, f)));
      }
    } else {
      T = state.temporal_size;
      const int c = T / 2;
      frames_holder.reserve(T);
      for (int j = 0; j < T; ++j) {
        const int real = (j >= c) ? ((N - 1 - n < j - c) ? (N - 1) : (n + (j - c)))
                                   : ((n < c - j) ? 0 : (n - (c - j)));
        frames_holder.push_back(unwrap(ctx.frames.get(0, real)));
      }
    }

    for (const auto& holder : frames_holder) {
      require(holder.frame.format == state.source.format && ctx.dst.format == state.source.format &&
                  holder.frame.plane_count == state.source.format.plane_count &&
                  ctx.dst.plane_count == holder.frame.plane_count,
              "frame format differs from plan");
    }

    if constexpr (A == Algorithm::FFT3D) {
      if (state.sampled) {
        for (int p = 0; p < state.source.format.plane_count; ++p) {
          const auto& plan = state.plans[p];
          if (!plan || plan->pattern_ready()) continue;
          const auto sample = unwrap(ctx.frames.get(0, state.pattern_frame));
          require(sample.frame.format == state.source.format && sample.frame.plane_count == state.source.format.plane_count,
                  "pattern frame format differs from plan");
          const auto& view = sample.frame.plane(p);
          switch (state.source.format.sample_format) {
            case ds::SampleFormat::UInt8: prepare_pattern<std::uint8_t>(view, *plan); break;
            case ds::SampleFormat::Float32: prepare_pattern<float>(view, *plan); break;
            default: prepare_pattern<std::uint16_t>(view, *plan); break;
          }
        }
      }
    }
    for (int p = 0; p < ctx.dst.plane_count; ++p) {
      const auto& d = ctx.dst.plane(p);
      const bool chroma = state.source.format.color_family == ds::ColorFamily::Yuv && p > 0;
      const int w = state.source.width >> (chroma ? state.source.format.subsampling_w : 0);
      const int h = state.source.height >> (chroma ? state.source.format.subsampling_h : 0);
      require(d.width == w && d.height == h, "frame plane dimensions differ from plan");

      std::vector<ds::PlaneView> plane_views;
      plane_views.reserve(T);
      for (int j = 0; j < T; ++j) {
        const auto& s = frames_holder[j].frame.plane(p);
        require(s.width == w && s.height == h, "frame plane dimensions differ from plan");
        plane_views.push_back(s);
      }

      switch (state.source.format.sample_format) {
        case ds::SampleFormat::UInt8:
          plane<std::uint8_t>(span2d::Span<const ds::PlaneView>(plane_views.data(), plane_views.size()), d,
                              state.plans[p].get());
          break;
        case ds::SampleFormat::Float32:
          plane<float>(span2d::Span<const ds::PlaneView>(plane_views.data(), plane_views.size()), d,
                       state.plans[p].get());
          break;
        default:
          plane<std::uint16_t>(span2d::Span<const ds::PlaneView>(plane_views.data(), plane_views.size()), d,
                               state.plans[p].get());
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
