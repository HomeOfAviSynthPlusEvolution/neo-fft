#pragma once
#include "plugin/descriptors.hpp"
#include "plugin/roi.hpp"
#include <dualsynth/staged_video.hpp>
#include "runtime/checkpoints.hpp"
#include <atomic>
#include <memory>
#include <cstring>

namespace neo_fft::plugin {
template <Algorithm A>
struct Filter {
  static constexpr const char* name = A == Algorithm::FFT3D ? "FFT3D" : "DFTTest";
  static constexpr int input_count = ds::dynamic_video_inputs;
  static constexpr ds::OutputOrigin output_origin = ds::OutputOrigin::fresh(0);
  static constexpr ds::HostRequirements host_requirements{true,0,0};
  struct State {
    ds::VideoInputInfo source;
    std::array<std::shared_ptr<const Plan>, 4> plans{};
    std::array<ROI, 4> rois{};
    std::array<PlaneMode, 4> plane_modes{PlaneMode::Copy, PlaneMode::Copy, PlaneMode::Copy, PlaneMode::Copy};
    int temporal_size = 1;
    int temporal_mode = 0, temporal_overlap = 0;
    int pattern_frame = 0;
    bool sampled = false;
    std::shared_ptr<const DFTNoise> dft_noise;
    int sample_bits = 8;
    bool kalman = false;
    std::shared_ptr<runtime::Retention> retention;
    std::shared_ptr<runtime::SpectraCache> spectra;
    std::shared_ptr<runtime::Checkpoints> checkpoints = std::make_shared<runtime::Checkpoints>();
    std::shared_ptr<std::atomic<bool>> models_ready = std::make_shared<std::atomic<bool>>(false);
  };
  static ds::Result<ds::VideoInitStateResult<State>> init(ds::VideoInitContext& ctx) {
    require(ctx.params && ctx.inputs.size() == 1, "missing clip or parameters");
    const auto& info = ctx.inputs[0];
    const auto& f = info.format;
    auto supported = ds::is_supported_video_format(f);
    require(supported.has_value() && supported.value() && (f.plane_count == 1 || f.plane_count == 3 || f.plane_count == 4),
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
    if constexpr (A == Algorithm::DFTTest) {
      require(config.tbsize <= info.num_frames, "tbsize must be less than or equal to the number of frames");
      for (const auto& location : config.locations) {
        require(location.frame <= info.num_frames-config.tbsize, "DFTTest sample interval outside clip");
        require(location.plane < f.plane_count, "DFTTest sample plane outside format");
        const bool chroma = f.color_family == ds::ColorFamily::Yuv && (location.plane == 1 || location.plane == 2);
        const int w=info.width >> (chroma ? f.subsampling_w : 0), h=info.height >> (chroma ? f.subsampling_h : 0);
        require(location.x <= w-config.block && location.y <= h-config.block, "DFTTest sample rectangle outside plane");
      }
    }
    const auto selected = select_planes(params, A, f.plane_count);
    const int t_size = [&] {
      if constexpr (A == Algorithm::FFT3D)
        return std::max(1, config.bt);
      else
        return config.tbsize;
    }();
    State state;
    state.plane_modes = selected;
    state.source = info; state.temporal_size = t_size;
    state.retention=std::make_shared<runtime::Retention>(1);
    if constexpr(A==Algorithm::DFTTest) {state.temporal_mode=config.temporal_mode;state.temporal_overlap=config.temporal_overlap;}
    state.sample_bits = ds::bits_per_sample(f.sample_format);
    for (int p = 0; p < f.plane_count; ++p)
      if (selected[p] == PlaneMode::Process) {
        try {
          const bool chroma = f.color_family == ds::ColorFamily::Yuv && (p == 1 || p == 2);
          const int w = info.width >> (chroma ? f.subsampling_w : 0), h = info.height >> (chroma ? f.subsampling_h : 0);
          const SampleFormat sample_format{state.sample_bits,f.sample_format == ds::SampleFormat::Float32,chroma};
          if constexpr (A == Algorithm::DFTTest) {
            state.plans[p] = std::make_shared<Plan>(w,h,sample_format,config,state.dft_noise,state.retention);
            state.dft_noise = state.plans[p]->dft_noise();
          } else {
            state.rois[p] = make_roi(w,h,chroma ? f.subsampling_w : 0,chroma ? f.subsampling_h : 0,config);
            state.plans[p] = std::make_shared<Plan>(state.rois[p].width,state.rois[p].height,sample_format,config,state.retention);
          }
        } catch (const std::exception& e) {
          throw std::invalid_argument("plane " + std::to_string(p) + ": " + e.what());
        }
      }
    if constexpr (A == Algorithm::FFT3D) {
      state.pattern_frame = std::clamp(config.pframe, 0, info.num_frames-1);
      std::array<std::size_t,4> bins{};
      std::array<int,4> rows{};
      for(int p=0;p<f.plane_count;++p)if(const auto& plan=state.plans[p]) {
        bins[p]=mul_size(plan->geometry.x.count,plan->fft.bins());
        rows[p]=plan->geometry.y.count;
      }
      if(config.bt>1)state.spectra=std::make_shared<runtime::SpectraCache>(bins,config.bt,config.cache_frames,config.cache_mb,rows);
      for (const auto& plan : state.plans) if (plan) {
        if (plan->preview()) state.temporal_size = 1;
        state.sampled = state.sampled || plan->needs_pattern_frame();
        state.kalman = state.kalman || plan->kalman();
      }
    }
    return ds::Result<ds::VideoInitStateResult<State>>::success(
        {{info.width, info.height, info.num_frames, f, info.fps}, std::move(state)});
  }
  static ds::VideoRequestPattern request_pattern(int, const State& state) {
    return state.kalman || state.temporal_size > 1 || state.sampled || state.dft_noise ? ds::VideoRequestPattern::General : ds::VideoRequestPattern::StrictSpatial;
  }
  static int cache_hints(ds::VideoCacheHintsContext& ctx) {
    // AviSynth V12 CACHE_INFORM_NUM_THREADS / AVS_CACHE_INFORM_NUM_THREADS.
    // Numeric ABI constant keeps the host-independent core free of AVS headers.
    if constexpr(A==Algorithm::FFT3D)if(ctx.cachehints==514) {
      auto& state=ctx.state<State>();
      if(state.spectra)state.spectra->inform_threads(ctx.frame_range);
    }
    return ctx.default_response;
  }
  struct RequestState {
    bool initialized=false, sample=false, ordinary=false;
    int pending=-1;
    std::unique_ptr<runtime::SpectraCache::Request> spectra;
    runtime::Checkpoints::Lease start;
    std::unique_ptr<runtime::Checkpoint> current;
    int replay_start=0, steps=0;
  };
  static std::unique_ptr<runtime::SpectraCache::Request> register_spectra(const State& state,int n) {
    if constexpr(A==Algorithm::FFT3D)if(state.spectra && state.temporal_size>1) {
      const int bt=state.temporal_size,left=bt/2,right=(bt-1)/2;
      if(n<left || state.source.num_frames-1-n<right)return state.spectra->register_request(n,n);
      return state.spectra->register_request(n-left,n+right);
    }
    return {};
  }
  template<class T,class Function> static void with_roi(const ds::PlaneView& view,const ROI& roi,CopyRow copy,Function&& fn) {
    auto s=checked_plane(static_cast<const T*>(view.data),view.width,view.height,view.stride_bytes,
                         plane_extent<T>(view.width,view.height,view.stride_bytes));
    if(!roi.interlaced) fn(checked_subplane(s,roi.left,roi.top,roi.width,roi.height));
    else { PackedROI<T> packed(s,roi,copy);fn(span2d::Plane<const T>(packed.view)); }
  }
  static void validate_source(const ds::VideoFrameView& frame,const State& state) {
    require(frame.format==state.source.format && frame.plane_count==state.source.format.plane_count,"frame format differs from plan");
    for(int p=0;p<frame.plane_count;++p) {
      const bool chroma=state.source.format.color_family==ds::ColorFamily::Yuv && (p==1 || p==2);
      require(frame.plane(p).width==(state.source.width>>(chroma ? state.source.format.subsampling_w:0)) &&
              frame.plane(p).height==(state.source.height>>(chroma ? state.source.format.subsampling_h:0)),"frame plane dimensions differ from plan");
    }
  }
  static ds::Result<ds::VideoStageResult> advance(ds::VideoStageContext& ctx,RequestState& r) {
    using Stage=ds::VideoStageResult;
    const auto& state=ctx.state<State>();const int n=ctx.output_frame;
    const auto waiting=[&](int index) {r.pending=index;ctx.request_frame(0,index);return ds::Result<Stage>::success(Stage::RequestFrames);};
    if(!r.initialized) {
      r.initialized=true;
      if(!state.kalman || n==0) {
        r.ordinary=true;
        if(state.kalman) ctx.request_frame(0,n);
        else {r.spectra=register_spectra(state,n);unwrap(request(ctx));}
        return ds::Result<Stage>::success(Stage::RequestFrames);
      }
      r.start=state.checkpoints->acquire(n);
      r.replay_start=r.start ? r.start->frame : 0;
      if(r.replay_start<n) {
        r.current=r.start ? std::make_unique<runtime::Checkpoint>(*r.start) : std::make_unique<runtime::Checkpoint>();
        if(!r.start) for(int p=0;p<state.source.format.plane_count;++p)
          if(state.plans[p]) r.current->planes[p]=state.plans[p]->initial_kalman();
      }
      r.sample=state.sampled && !state.models_ready->load(std::memory_order_acquire);
      if(r.sample) return waiting(state.pattern_frame);
      return waiting(r.replay_start<n ? r.replay_start+1 : n);
    }
    if(r.ordinary) return ds::Result<Stage>::success(Stage::Ready);
    {
      const auto source=unwrap(ctx.frames.get(0,r.pending));
      validate_source(source.frame,state);
      if(r.sample) {
        std::array<runtime::PublishedModel::Model,4> candidates;
        for(int p=0;p<state.source.format.plane_count;++p) if(state.plans[p] && state.plans[p]->needs_pattern_frame()) {
          const auto compute=[&](auto s) {candidates[p]=state.plans[p]->pattern_candidate(s);};
          const auto& view=source.frame.plane(p);
          switch(state.source.format.sample_format) {
            case ds::SampleFormat::UInt8: with_roi<std::uint8_t>(view,state.rois[p],state.plans[p]->copy_row(),compute);break;
            case ds::SampleFormat::Float32: with_roi<float>(view,state.rois[p],state.plans[p]->copy_row(),compute);break;
            default: with_roi<std::uint16_t>(view,state.rois[p],state.plans[p]->copy_row(),compute);break;
          }
        }
        for(int p=0;p<state.source.format.plane_count;++p) if(candidates[p]) state.plans[p]->publish_pattern(std::move(candidates[p]));
        state.models_ready->store(true,std::memory_order_release);
      } else if(r.current) {
        for(int p=0;p<state.source.format.plane_count;++p) { if(!state.plans[p]) continue;
          const auto consume=[&](auto s) {state.plans[p]->advance_kalman(s,r.current->planes[p]);};
          const auto& view=source.frame.plane(p);
          switch(state.source.format.sample_format) {
            case ds::SampleFormat::UInt8: with_roi<std::uint8_t>(view,state.rois[p],state.plans[p]->copy_row(),consume);break;
            case ds::SampleFormat::Float32: with_roi<float>(view,state.rois[p],state.plans[p]->copy_row(),consume);break;
            default: with_roi<std::uint16_t>(view,state.rois[p],state.plans[p]->copy_row(),consume);break;
          }
        }
        r.current->frame=r.pending; ++r.steps;
      }
    } // Drop get()'s owning snapshot before releasing the staged-store owner.
    if(r.sample) {
      unwrap(ctx.release_frame(0,r.pending));r.sample=false;
      return waiting(r.replay_start<n ? r.replay_start+1 : n);
    }
    if(r.pending==n) return ds::Result<Stage>::success(Stage::Ready);
    unwrap(ctx.release_frame(0,r.pending));
    return waiting(r.pending+1); // pending<n<=INT_MAX, so this cannot overflow.
  }
  template<class T> static void render_plane(const ds::PlaneView& src,const ds::MutablePlaneView& dst,
                                             const Plan* plan,const ROI& roi,const KalmanState* state) {
    const auto s=checked_plane(static_cast<const T*>(src.data),src.width,src.height,src.stride_bytes,plane_extent<T>(src.width,src.height,src.stride_bytes));
    const auto d=checked_plane(static_cast<T*>(dst.data),dst.width,dst.height,dst.stride_bytes,plane_extent<T>(dst.width,dst.height,dst.stride_bytes));
    require(s.width()==d.width() && s.height()==d.height(),"output plane shape differs");
    disjoint(s.data(),plane_extent<T>(src.width,src.height,src.stride_bytes),d.data(),plane_extent<T>(dst.width,dst.height,dst.stride_bytes));
    for(int y=0;y<d.height();++y) std::memcpy(d.row_ptr(y),s.row_ptr(y),std::size_t(d.width())*sizeof(T));
    if(!plan || !state) return;
    if(!roi.interlaced) plan->render_kalman(checked_subplane(s,roi.left,roi.top,roi.width,roi.height),checked_subplane(d,roi.left,roi.top,roi.width,roi.height),*state);
    else {
      PackedROI<T> input(s,roi,plan->copy_row()),output(s,roi,plan->copy_row());
      plan->render_kalman(span2d::Plane<const T>(input.view),output.view,*state);output.write(d,roi);
    }
  }
  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& ctx,RequestState& r) {
    const auto& state=ctx.state<State>();
    if(!state.kalman) {
      auto registration=std::move(r.spectra);
      return process(ctx,registration.get());
    }
    const auto src=unwrap(ctx.frames.get(0,ctx.output_frame));validate_source(src.frame,state);
    require(ctx.dst.format==state.source.format && ctx.dst.plane_count==state.source.format.plane_count,"output format differs");
    const auto* checkpoint=r.current ? r.current.get() : r.start.get();
    for(int p=0;p<ctx.dst.plane_count;++p) {
      if(state.plane_modes[p]==PlaneMode::Skip) continue;
      const auto* k=checkpoint && state.plans[p] ? &checkpoint->planes[p] : nullptr;
      const auto& s=src.frame.plane(p);const auto& d=ctx.dst.plane(p);
      switch(state.source.format.sample_format) {
        case ds::SampleFormat::UInt8: render_plane<std::uint8_t>(s,d,state.plans[p].get(),state.rois[p],k);break;
        case ds::SampleFormat::Float32: render_plane<float>(s,d,state.plans[p].get(),state.rois[p],k);break;
        default: render_plane<std::uint16_t>(s,d,state.plans[p].get(),state.rois[p],k);break;
      }
    }
    if(r.current) state.checkpoints->publish(std::move(r.current));
    return ds::Result<ds::VideoProcessResult>::success({});
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
      std::vector<int> requested;
      if (state.dft_noise) for (const auto& location : state.dft_noise->locations)
        for (int j=0;j<state.temporal_size;++j) requested.push_back(location.frame+j);
      const int T = state.temporal_size;
      if(state.temporal_mode==1) {
        for(const auto& block:temporal_blocks(n,N,T,state.temporal_overlap))
          for(int z=0;z<T;++z)requested.push_back(block.slots[z]);
      } else {
        const int c = T / 2;
        for (int j = 0; j < T; ++j) {
          const int real = (j >= c) ? ((N - 1 - n < j - c) ? (N - 1) : (n + (j - c)))
                                     : ((n < c - j) ? 0 : (n - (c - j)));
          requested.push_back(real);
        }
      }
      std::sort(requested.begin(),requested.end());
      requested.erase(std::unique(requested.begin(),requested.end()),requested.end());
      if(ctx.requests.empty()) {
        // DualSynth's request_frame performs a linear duplicate search. The
        // sorted set is already unique, so append it directly in the common path.
        ctx.requests.reserve(requested.size());
        for(int frame:requested)ctx.requests.push_back({0,frame});
      } else {
        for(int frame:requested)ctx.request_frame(0,frame);
      }
    }
    return ds::Result<ds::VideoRequestResult>::success({});
  }
  template <class T>
  static void plane(span2d::Span<const ds::PlaneView> src_views, const ds::MutablePlaneView& dst, const Plan* plan, const ROI& roi, int frame, int plane_index,span2d::Span<const int> targets={},runtime::SpectraCache* spectra_cache=nullptr,runtime::SpectraCache::Request* registration=nullptr) {
    const auto de = plane_extent<T>(dst.width, dst.height, dst.stride_bytes);
    auto d = checked_plane(static_cast<T*>(static_cast<void*>(dst.data)), dst.width, dst.height, dst.stride_bytes, de);

    if (plan) {
      require(src_views.size() <= 225, "too many temporal slots");
      std::vector<span2d::Plane<const T>> checked_srcs;
      checked_srcs.reserve(src_views.size());
      // A fixed open-addressed index bounds lookup even at T=15/O=14,
      // where 225 logical slots refer to at most 29 output-neighborhood frames.
      std::array<std::size_t,512> cache{}; // 0 is empty; otherwise slot index + 1.
      for (std::size_t i=0;i<src_views.size();++i) {
        const auto& src=src_views[i];
        std::size_t slot=(reinterpret_cast<std::uintptr_t>(src.data)>>4)&(cache.size()-1);
        bool seen=false;
        while(cache[slot]) {
          const auto prior_index=cache[slot]-1;
          const auto& prior=src_views[prior_index];
          if(src.data==prior.data && src.width==prior.width && src.height==prior.height &&
             src.stride_bytes==prior.stride_bytes) {
            checked_srcs.push_back(checked_srcs[prior_index]);
            seen=true;
            break;
          }
          slot=(slot+1)&(cache.size()-1);
        }
        if(seen) continue;
        const auto se=plane_extent<T>(src.width,src.height,src.stride_bytes);
        checked_srcs.push_back(checked_plane(static_cast<const T*>(static_cast<const void*>(src.data)),
                                             src.width,src.height,src.stride_bytes,se));
        disjoint(checked_srcs.back().data(),se,d.data(),de);
        cache[slot]=i+1;
      }
      if constexpr (A == Algorithm::FFT3D) {
        const auto& original = checked_srcs[checked_srcs.size()/2];
        for (int y=0;y<dst.height;++y)
          std::memcpy(d.row_ptr(y),original.row_ptr(y),std::size_t(dst.width)*sizeof(T));
        if (!roi.interlaced) {
          for (auto& source : checked_srcs)
            source = checked_subplane(source,roi.left,roi.top,roi.width,roi.height);
          plan->process_at<T>({checked_srcs.data(),checked_srcs.size()},checked_subplane(d,roi.left,roi.top,roi.width,roi.height),frame,plane_index,{},spectra_cache,registration);
          return;
        }
        std::vector<PackedROI<T>> packed;
        packed.reserve(checked_srcs.size());
        for (auto source : checked_srcs) packed.emplace_back(source,roi,plan->copy_row());
        std::vector<span2d::Plane<const T>> views;
        for (auto& image : packed) views.push_back(image.view);
        PackedROI<T> output(original,roi,plan->copy_row());
        plan->process_at<T>({views.data(),views.size()},output.view,frame,plane_index,{},spectra_cache,registration);
        output.write(d,roi);
      } else plan->process_at<T>({checked_srcs.data(),checked_srcs.size()},d,frame,plane_index,targets);
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
  template<class T> static void prepare_pattern(const ds::PlaneView& src, const Plan& plan, const ROI& roi) {
    const auto s = checked_plane(static_cast<const T*>(static_cast<const void*>(src.data)), src.width, src.height,
        src.stride_bytes, plane_extent<T>(src.width, src.height, src.stride_bytes));
    if (!roi.interlaced) plan.prepare_pattern(checked_subplane(s,roi.left,roi.top,roi.width,roi.height));
    else {
      PackedROI<T> packed(s,roi,plan.copy_row());
      plan.prepare_pattern(packed.view);
    }
  }
  template<class T> static void gather_noise(const ds::PlaneView& view, const NoiseLocation& location,
                                            int S, int bits, span2d::Span<float> out, const ModelKernels& kernels) {
    const auto source = checked_plane(static_cast<const T*>(static_cast<const void*>(view.data)),view.width,view.height,
        view.stride_bytes,plane_extent<T>(view.width,view.height,view.stride_bytes));
    require(location.x <= view.width-S && location.y <= view.height-S, "sample frame rectangle differs from plan");
    const float scale=std::is_same_v<T,float> ? 255.0f : 1.0f/float(1 << (bits-8));
    for(int y=0;y<S;++y)
      kernels.decode(source.row_ptr(location.y+y)+location.x,sample_storage<T>,out.data()+std::size_t(y)*S,std::size_t(S),0,scale);
  }
  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& ctx,runtime::SpectraCache::Request* registration=nullptr) {
    const auto& state = ctx.state<State>();
    const int n = ctx.output_frame;
    const int N = state.source.num_frames;
    auto local_registration=registration ? nullptr : register_spectra(state,n);
    if(!registration)registration=local_registration.get();

    const bool has_active_plans = std::any_of(state.plans.begin(), state.plans.end(), [](const auto& p) {
      return p != nullptr;
    });

    int T = 1;
    std::vector<ds::RequestedVideoFrame> frames_holder;
    std::vector<int> targets;
    std::vector<std::size_t> slot_to_frame;

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
      std::vector<int> slot_frames;
      if(state.temporal_mode==1) {
        for(const auto& block:temporal_blocks(n,N,T,state.temporal_overlap)) {
          targets.push_back(block.target);
          for(int z=0;z<T;++z)slot_frames.push_back(block.slots[z]);
        }
      } else {
      const int c = T / 2;
      slot_frames.reserve(T);
      for (int j = 0; j < T; ++j) {
        const int real = (j >= c) ? ((N - 1 - n < j - c) ? (N - 1) : (n + (j - c)))
                                   : ((n < c - j) ? 0 : (n - (c - j)));
        slot_frames.push_back(real);
      }
      }
      auto unique_frames=slot_frames;
      std::sort(unique_frames.begin(),unique_frames.end());
      unique_frames.erase(std::unique(unique_frames.begin(),unique_frames.end()),unique_frames.end());
      frames_holder.reserve(unique_frames.size());
      for(int frame:unique_frames)frames_holder.push_back(unwrap(ctx.frames.get(0,frame)));
      slot_to_frame.reserve(slot_frames.size());
      for(int frame:slot_frames)slot_to_frame.push_back(std::size_t(std::lower_bound(unique_frames.begin(),unique_frames.end(),frame)-unique_frames.begin()));
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
          const bool chroma = state.source.format.color_family == ds::ColorFamily::Yuv && (p == 1 || p == 2);
          require(view.width == (state.source.width >> (chroma ? state.source.format.subsampling_w : 0)) &&
                  view.height == (state.source.height >> (chroma ? state.source.format.subsampling_h : 0)),
                  "pattern plane dimensions differ from plan");
          switch (state.source.format.sample_format) {
            case ds::SampleFormat::UInt8: prepare_pattern<std::uint8_t>(view, *plan, state.rois[p]); break;
            case ds::SampleFormat::Float32: prepare_pattern<float>(view, *plan, state.rois[p]); break;
            default: prepare_pattern<std::uint16_t>(view, *plan, state.rois[p]); break;
          }
        }
      }
    }
    if constexpr (A == Algorithm::DFTTest) {
      if (state.dft_noise) state.dft_noise->prepare([&](const NoiseLocation& location,int z,span2d::Span<float> out) {
        const auto sample=unwrap(ctx.frames.get(0,location.frame+z));
        require(sample.frame.format==state.source.format && sample.frame.plane_count==state.source.format.plane_count,
                "sample frame format differs from plan");
        const auto& view=sample.frame.plane(location.plane);
        const bool chroma=state.source.format.color_family==ds::ColorFamily::Yuv && (location.plane==1 || location.plane==2);
        require(view.width==(state.source.width >> (chroma ? state.source.format.subsampling_w : 0)) &&
                view.height==(state.source.height >> (chroma ? state.source.format.subsampling_h : 0)), "sample plane dimensions differ from plan");
        const int S=state.dft_noise->block_size, bits=state.sample_bits;
        switch (state.source.format.sample_format) {
          case ds::SampleFormat::UInt8: gather_noise<std::uint8_t>(view,location,S,bits,out,state.dft_noise->kernels()); break;
          case ds::SampleFormat::Float32: gather_noise<float>(view,location,S,bits,out,state.dft_noise->kernels()); break;
          default: gather_noise<std::uint16_t>(view,location,S,bits,out,state.dft_noise->kernels()); break;
        }
      });
    }
    const auto process_plane = [&](int p) {
      if(state.plane_modes[p]==PlaneMode::Skip) return;
      const auto& d = ctx.dst.plane(p);
      const bool chroma = state.source.format.color_family == ds::ColorFamily::Yuv && (p == 1 || p == 2);
      const int w = state.source.width >> (chroma ? state.source.format.subsampling_w : 0);
      const int h = state.source.height >> (chroma ? state.source.format.subsampling_h : 0);
      require(d.width == w && d.height == h, "frame plane dimensions differ from plan");

      std::vector<ds::PlaneView> plane_views;
      plane_views.reserve(slot_to_frame.empty() ? frames_holder.size() : slot_to_frame.size());
      for (std::size_t j = 0; j < frames_holder.size(); ++j) {
        const auto& s = frames_holder[j].frame.plane(p);
        require(s.width == w && s.height == h, "frame plane dimensions differ from plan");
        if(slot_to_frame.empty())plane_views.push_back(s);
      }
      for(std::size_t index:slot_to_frame)plane_views.push_back(frames_holder[index].frame.plane(p));

      // The midpoint of concatenated block slots need not be frame n.
      if(!state.plans[p] && !targets.empty())plane_views={frames_holder[slot_to_frame[std::size_t(targets[0])]].frame.plane(p)};
      switch (state.source.format.sample_format) {
        case ds::SampleFormat::UInt8:
          plane<std::uint8_t>(span2d::Span<const ds::PlaneView>(plane_views.data(), plane_views.size()), d,
                              state.plans[p].get(),state.rois[p],n,p,{targets.data(),targets.size()},state.spectra.get(),registration);
          break;
        case ds::SampleFormat::Float32:
          plane<float>(span2d::Span<const ds::PlaneView>(plane_views.data(), plane_views.size()), d,
                       state.plans[p].get(),state.rois[p],n,p,{targets.data(),targets.size()},state.spectra.get(),registration);
          break;
        default:
          plane<std::uint16_t>(span2d::Span<const ds::PlaneView>(plane_views.data(), plane_views.size()), d,
                               state.plans[p].get(),state.rois[p],n,p,{targets.data(),targets.size()},state.spectra.get(),registration);
          break;
      }
    };
    for(int p=0;p<ctx.dst.plane_count;++p) process_plane(p);
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
    return v.has_value() && v.value() && (f.plane_count == 1 || f.plane_count == 3 || f.plane_count == 4);
  }
};
} // namespace neo_fft::plugin
