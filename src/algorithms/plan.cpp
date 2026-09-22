#include "algorithms/plan.hpp"
#include "algorithms/pad.hpp"
#include "kernels/spectral.hpp"
#include <cstring>
#include <type_traits>

namespace neo_fft {
namespace {
void nonnegative(float v, const char* name) {
  require(std::isfinite(v) && v >= 0, std::string("invalid ") + name);
}
void valid_format(SampleFormat f) {
  require(f.floating ? f.bits == 32 : (f.bits == 8 || f.bits == 10 || f.bits == 12 || f.bits == 14 || f.bits == 16),
          "unsupported sample format");
}
Geometry geometry3d(int w, int h, const FFT3DConfig& c) {
  validate(c);
  return {fft3d_axis(w, c.bw, c.ow), fft3d_axis(h, c.bh, c.oh)};
}
Geometry geometry_dft(int w, int h, const DFTConfig& c) {
  validate(c);
  return {dft_axis(w, c.block, c.mode, c.overlap), dft_axis(h, c.block, c.mode, c.overlap)};
}
int select_optimal_batch_size(std::size_t block_samples, int gx_count) noexcept {
  const std::size_t l2_budget = optimal_l2_working_set_bytes();
  const int lanes = std::max(4, optimal_simd_lanes());
  const std::size_t bytes_per_block = block_samples * 16;
  int k_max = static_cast<int>(l2_budget / (bytes_per_block ? bytes_per_block : 1));
  k_max = (k_max / lanes) * lanes;
  k_max = std::clamp(k_max, lanes, 64);

  int best_k = k_max;
  for (int cand = k_max; cand >= std::max(lanes, k_max / 2); cand -= lanes) {
    if (gx_count % cand == 0) {
      best_k = cand;
      break;
    }
  }
  return std::max(1, std::min(gx_count, best_k));
}
} // namespace
void validate(const FFT3DConfig& c) {
  require(c.bw >= 2 && c.bh >= 2, "FFT3D bw/bh must be >=2");
  require(c.ow <= c.bw / 2 && c.oh <= c.bh / 2, "FFT3D ow/oh exceeds half block");
  require(c.wintype >= 0 && c.wintype <= 2, "FFT3D wintype outside 0..2");
  require(c.bt == -1 || (c.bt >= 1 && c.bt <= 5), "FFT3D bt outside -1 or 1..5");
  validate(c.enhancement);
  nonnegative(c.pfactor, "FFT3D pfactor");
  require(std::isfinite(c.pcutoff) && c.pcutoff > 0, "FFT3D pcutoff must be positive");
  require(c.px >= 0 && c.py >= 0, "FFT3D px/py must be nonnegative");
  for (auto s : {c.sigma2, c.sigma3, c.sigma4}) nonnegative(s.value_or(c.sigma), "FFT3D sigma2..4");
  nonnegative(c.sigma, "FFT3D sigma");
  nonnegative(c.degrid, "FFT3D degrid");
  require(std::isfinite(c.beta) && c.beta >= 1, "FFT3D beta must be >=1");
  select_spectral(c.opt);
}
void validate(const DFTConfig& c) {
  validate(c.curves);
  require(c.locations.size() <= 500, "DFTTest nlocation exceeds 500 tuples");
  const float alpha = c.alpha.value_or(c.ftype == 0 ? 5.0f : 7.0f);
  require(std::isfinite(alpha) && alpha > 0, "DFTTest alpha must be positive");
  for (const auto& n : c.locations) require(n.frame >= 0 && n.plane >= 0 && n.y >= 0 && n.x >= 0, "DFTTest negative sample coordinate");
  require(c.block > 0 && (c.mode == 0 || c.mode == 1), "DFTTest invalid sbsize/smode");
  require(c.mode != 0 || c.block % 2 == 1, "DFTTest center mode requires odd sbsize");
  require(c.mode == 0 || (c.overlap >= 0 && c.overlap < c.block), "DFTTest invalid sosize");
  require(c.tbsize >= 1 && c.tbsize <= 15 && c.tbsize % 2 == 1, "DFTTest tbsize must be odd integer in 1..15");
  if (c.mode == 1 && c.overlap > c.block / 2)
    require(c.block % (c.block - c.overlap) == 0, "DFTTest heavy overlap requires divisible step");
  require(c.swin >= 0 && c.swin <= 11 && c.twin >= 0 && c.twin <= 11, "DFTTest invalid swin/twin");
  require(c.ftype >= 0 && c.ftype <= 4, "DFTTest ftype outside 0..4");
  nonnegative(c.sbeta, "DFTTest sbeta");
  nonnegative(c.tbeta, "DFTTest tbeta");
  nonnegative(c.sigma, "DFTTest sigma");
  nonnegative(c.sigma2, "DFTTest sigma2");
  nonnegative(c.pmin, "DFTTest pmin");
  nonnegative(c.pmax, "DFTTest pmax");
  require(c.pmin <= c.pmax, "DFTTest pmin exceeds pmax");
  require(std::isfinite(c.f0beta) && c.f0beta > 0, "DFTTest f0beta must be positive");
  select_spectral(c.opt);
  select_spatial(c.opt);
}
Plan::Plan(int w, int h, SampleFormat f, const FFT3DConfig& c)
    : geometry(geometry3d(w, h, c)), format(f), algorithm(Algorithm::FFT3D), temporal_size(std::max(1, c.bt)), fft(c.bh, c.bw),
      denoise_(c.bt != -1),
      wx_(fft3d_window(c.bw, geometry.x.overlap, c.wintype)), wy_(fft3d_window(c.bh, geometry.y.overlap, c.wintype)),
      kernel_(select_spectral(c.opt)), temporal_kernel_(select_fft3d_temporal(c.opt)),
      spatial_(select_spatial(c.opt)), mean_scale_(c.degrid),
      pool_(runtime::make_workspace_budget(geometry, fft.samples(),
                                           fft.bins() * std::size_t(std::max(1, c.bt) + 1), true,
                                           std::max(1, c.bt),
                                           select_optimal_batch_size(fft.samples(), geometry.x.count))) {
  valid_format(f);
  const float factor = f.floating ? 1.0f / 255 : float(1 << (f.bits - 8));
  sampled_ = c.pfactor > 0;
  pfactor_ = c.pfactor;
  px_ = c.px; py_ = c.py;
  sigma_eff_ = sampled_ ? 0 : finite(c.sigma * factor);
  norm_ = 1.0f / (float(c.bw) * float(c.bh));
  beta_ = c.beta;
  const std::array<float, 4> sigmas = sampled_ ? std::array<float,4>{} : std::array<float,4>{sigma_eff_, finite(c.sigma2.value_or(c.sigma) * factor),
      finite(c.sigma3.value_or(c.sigma) * factor), finite(c.sigma4.value_or(c.sigma) * factor)};
  const bool varying = sigmas[1] != sigmas[0] || sigmas[2] != sigmas[0] || sigmas[3] != sigmas[0];
  preview_ = c.pshow && (sampled_ || varying);
  denoise_ = denoise_ && !preview_;
  params_.floor = (c.beta - 1) / c.beta;
  if (denoise_ && !sampled_) {
    if (varying) {
      primary_ = fft3d_profile(c.bw, c.bh, sigmas);
      params_.primary_mode = PrimaryMode::Table;
      params_.primary = {primary_.data(), primary_.size()};
      for (float v : primary_) finite(float(temporal_size) * v);
    } else {
      params_.a = finite((sigma_eff_ * sigma_eff_) / norm_);
      for (int T = 1; T <= temporal_size; ++T) finite(((float(T) * sigma_eff_) * sigma_eff_) / norm_);
    }
  }
  if (!preview_) enhancement_tables_ = enhancement_tables(c.bw, c.bh, factor, c.enhancement);
  enhancement_params_.type = -2;
  if (!preview_) enhancement_params_.enhancement = enhancement_tables_.view(c.enhancement);
  if (!varying && !sampled_ && denoise_) params_.enhancement = enhancement_params_.enhancement;
  if (preview_ || needs_pattern_frame()) {
    if (px_ == 0 && py_ == 0) require(geometry.x.count >= 5 && geometry.y.count >= 5, "FFT3D automatic pattern search requires a 5x5 grid");
    else require(px_ < geometry.x.count && py_ < geometry.y.count, "FFT3D pattern coordinates outside grid");
    sample_weights_ = buffer<float>(fft.bins());
    const float cutoff2 = finite(c.pcutoff * c.pcutoff);
    require(cutoff2 > 0, "FFT3D pattern cutoff underflow");
    for (int y = 0; y < c.bh; ++y) for (int x = 0; x <= c.bw / 2; ++x) {
      const float fy = 2 * float(std::min(y, c.bh-y)) / float(c.bh), fx = 2 * float(x) / float(c.bw);
      const float q = fy*fy + fx*fx;
      sample_weights_[std::size_t(y)*(c.bw/2+1)+x] = finite(q / finite(q+cutoff2));
    }
  }
  if (c.degrid != 0) {
    auto block = buffer<float>(fft.samples());
    const float peak = f.floating ? 1.0f : float((1 << f.bits) - 1);
    for (int y = 0; y < c.bh; ++y)
      for (int x = 0; x < c.bw; ++x)
        block[std::size_t(y) * c.bw + x] = (peak * wy_.analysis[y]) * wx_.analysis[x];
    grid_ = buffer<std::complex<float>>(fft.bins());
    fft.forward(block.data(), grid_.data());
    require(grid_[0].real() != 0, "FFT3D unusable grid DC");
  }
}
Plan::Plan(int w, int h, SampleFormat f, const DFTConfig& c, std::shared_ptr<const DFTNoise> noise)
    : geometry(geometry_dft(w, h, c)), format(f), algorithm(Algorithm::DFTTest), temporal_size(c.tbsize),
      fft(c.block, c.block), kernel_(select_spectral(c.opt)), spatial_(select_spatial(c.opt)),
      mean_scale_(c.zmean ? 1.0f : 0.0f), center_(c.mode == 0),
      pool_(runtime::make_workspace_budget(
          geometry,
          std::size_t(c.tbsize) * c.block * c.block,
          std::size_t(c.tbsize) * c.block * (c.block / 2 + 1),
          false,
          c.tbsize,
          select_optimal_batch_size(std::size_t(c.tbsize) * c.block * c.block, geometry.x.count))) {
  valid_format(f);
  if (c.tbsize > 1) {
    fft3d_ = std::make_unique<RealFFT3D>(c.tbsize, c.block, c.block);
  }
  auto win = dft_window_3d(c.tbsize, c.block, c.mode == 0 ? 0 : c.overlap, c.mode, c.swin, c.twin, c.sbeta, c.tbeta);
  h_ = std::move(win.h);
  h_synthesis_.resize(h_.size());
  const float volume = float(c.tbsize) * float(c.block) * float(c.block);
  for (std::size_t i = 0; i < h_.size(); ++i)
    h_synthesis_[i] = h_[i] * volume;
  const bool sampled = c.ftype < 2 && !c.locations.empty();
  if (sampled) dft_noise_ = noise ? std::move(noise) : std::make_shared<DFTNoise>(c);
  const float scale = c.ftype < 2 ? win.wscale : 1.0f;
  params_ = {c.ftype,
             !sampled && c.curves.empty() ? finite(c.sigma / scale) : 0,
             finite(c.sigma2 / scale),
             finite(c.pmin / win.wscale),
             finite(c.pmax / win.wscale),
             c.f0beta,
             0};
  if (!sampled && !c.curves.empty()) {
    primary_ = dft_profile(c.curves, c.tbsize, c.block, c.sigma, scale);
    params_.primary_mode = PrimaryMode::Table;
    params_.primary = {primary_.data(), primary_.size()};
  }
  if (c.zmean) {
    auto block = h_;
    for (float& v : block)
      v = finite(255.0f * v);
    if (fft3d_) {
      grid_ = buffer<std::complex<float>>(fft3d_->bins());
      fft3d_->forward(block.data(), grid_.data());
    } else {
      grid_ = buffer<std::complex<float>>(fft.bins());
      fft.forward(block.data(), grid_.data());
    }
    require(grid_[0].real() != 0, "DFTTest zmean requires nonzero template DC");
  }
}
template<class T>
std::pair<int,int> Plan::select_pattern(span2d::Plane<const T> source, runtime::Workspace& ws, std::vector<float>* power) const {
  require(source.width() == geometry.x.length && source.height() == geometry.y.length, "pattern plane dimensions differ from plan");
  require(format.floating == std::is_same_v<T,float> && (format.floating || ((format.bits == 8) == (sizeof(T) == 1))), "pattern storage differs from plan");
  checked_plane(source.data(), source.width(), source.height(), source.stride_bytes(),
                plane_extent<T>(source.width(), source.height(), source.stride_bytes()));
  const bool automatic = px_ == 0 && py_ == 0;
  const int left = automatic ? 2 : px_, right = automatic ? geometry.x.count-3 : px_;
  const int top = automatic ? 2 : py_, bottom = automatic ? geometry.y.count-3 : py_;
  float best = 0;
  bool found = false;
  std::pair<int,int> selected{};
  const float base = !format.floating && format.chroma ? float(1 << (format.bits-1)) : 0;
  for (int by = top; by <= bottom; ++by) for (int bx = left; bx <= right; ++bx) {
    for (int y = 0; y < geometry.y.block; ++y) for (int x = 0; x < geometry.x.block; ++x) {
      const int sy = reflect(by*geometry.y.step+y, geometry.y), sx = reflect(bx*geometry.x.step+x, geometry.x);
      const float value = finite(finite(float(source.row_ptr(sy)[sx])) - base);
      ws.block()[std::size_t(y)*geometry.x.block+x] = finite(finite(value*wy_.analysis[y])*wx_.analysis[x]);
    }
    auto* spectrum = ws.spectrum().data();
    fft.forward(ws.block().data(), spectrum);
    const float g = grid_.empty() ? 0 : finite(finite(mean_scale_*spectrum[0].real())/grid_[0].real());
    float score = 0;
    for (std::size_t k = 0; k < fft.bins(); ++k) {
      const float re = finite(spectrum[k].real() - (grid_.empty() ? 0 : finite(g*grid_[k].real())));
      const float im = finite(spectrum[k].imag() - (grid_.empty() ? 0 : finite(g*grid_[k].imag())));
      const float q = finite(finite(re*re) + finite(im*im));
      // Request-private inverse storage is available for temporary unscaled powers.
      ws.inverse()[k] = q;
      score = finite(score + finite(q * sample_weights_[k]));
    }
    if (!found || score < best) {
      found = true; best = score; selected = {bx,by};
      if (power) for (std::size_t k = 0; k < fft.bins(); ++k) (*power)[k] = ws.inverse()[k];
    }
  }
  require(found, "FFT3D no pattern candidate");
  return selected;
}
template<class T> void Plan::build_pattern(span2d::Plane<const T> source) const {
  require(needs_pattern_frame(), "FFT3D sampling is inactive");
  if (sampled_model_.get()) return;
  auto lease = pool_.acquire();
  auto power = buffer<float>(fft.bins());
  select_pattern(source, *lease, &power);
  for (std::size_t k = 0; k < power.size(); ++k) {
    power[k] = finite(finite(pfactor_*power[k])*sample_weights_[k]);
    finite(float(temporal_size)*power[k]);
  }
  sampled_model_.publish(std::move(power));
}
void Plan::prepare_pattern(span2d::Plane<const std::uint8_t> source) const { build_pattern(source); }
void Plan::prepare_pattern(span2d::Plane<const std::uint16_t> source) const { build_pattern(source); }
void Plan::prepare_pattern(span2d::Plane<const float> source) const { build_pattern(source); }

void Plan::enhance(std::complex<float>* spectrum) const {
  if (!enhancement_params_.enhancement.active()) return;
  const float scale = grid_.empty() ? 0 : finite(finite(mean_scale_ * spectrum[0].real()) / grid_[0].real());
  kernel_(spectrum, grid_.empty() ? nullptr : grid_.data(), fft.bins(), scale, enhancement_params_);
}

template <class T>
void Plan::run(span2d::Span<const span2d::Plane<const T>> sources, span2d::Plane<T> dst, runtime::Workspace& ws) const {
  require(!sources.empty(), "sources must not be empty");
  const auto& gx = geometry.x;
  const auto& gy = geometry.y;
  const int T_slots = int(sources.size());
  if (algorithm == Algorithm::DFTTest) {
    require(T_slots == temporal_size, "DFTTest sources size must match plan temporal_size");
  } else {
    require(T_slots == temporal_size || T_slots == 1,
            "FFT3D sources size must match plan temporal_size or single-frame fallback");
  }
  require(ws.budget().temporal_slots >= T_slots, "workspace temporal slots insufficient");
  for (int j = 0; j < T_slots; ++j) {
    const auto& src = sources[j];
    require(src.width() == gx.length && src.height() == gy.length && dst.width() == gx.length &&
                dst.height() == gy.length,
            "frame dimensions differ from plan");
    require(format.floating == std::is_same_v<T, float> && (format.floating || ((format.bits == 8) == (sizeof(T) == 1))),
            "frame storage differs from plan");
    const auto se = plane_extent<T>(src.width(), src.height(), src.stride_bytes());
    checked_plane(src.data(), src.width(), src.height(), src.stride_bytes(), se);
    if constexpr (std::is_same_v<T, float>) {
      for (int y = 0; !preview_ && y < src.height(); ++y)
        for (int x = 0; x < src.width(); ++x)
          finite(src.row_ptr(y)[x]);
    }
  }
  const auto de = plane_extent<T>(dst.width(), dst.height(), dst.stride_bytes());
  checked_plane(dst.data(), dst.width(), dst.height(), dst.stride_bytes(), de);
  for (int j = 0; j < T_slots; ++j) {
    const auto se = plane_extent<T>(sources[j].width(), sources[j].height(), sources[j].stride_bytes());
    disjoint(sources[j].data(), se, dst.data(), de);
  }

  auto parameters = params_;
  const auto sampled_model = needs_pattern_frame() ? sampled_model_.get() : runtime::PublishedModel::Model{};
  if (needs_pattern_frame()) {
    require(bool(sampled_model), "FFT3D sampled model not initialized");
    parameters.primary_mode = PrimaryMode::Table;
    parameters.primary = {sampled_model->data(), sampled_model->size()};
  }
  const auto dft_model = dft_noise_ ? dft_noise_->power() : runtime::PublishedModel::Model{};
  if (dft_noise_) {
    require(bool(dft_model), "DFTTest sample model not initialized");
    parameters.primary_mode = PrimaryMode::Table;
    parameters.primary = {dft_model->data(),dft_model->size()};
  }
  ws.reset();
  for (int j = 0; j < T_slots; ++j) {
    if (!preview_) pad_source(sources[j], ws.padded(j), geometry, format, algorithm);
  }

  auto accum = ws.accum();
  auto row = algorithm == Algorithm::FFT3D ? ws.row() : span2d::Plane<float>{};
  const float base =
      !format.floating && format.chroma && algorithm == Algorithm::FFT3D ? float(1 << (format.bits - 1)) : 0;

  if (algorithm == Algorithm::FFT3D && preview_) {
    const auto selected = select_pattern(sources[0], ws, nullptr);
    const int ox = selected.first * gx.step, oy = selected.second * gy.step;
    for (int y = 0; y < gy.block; ++y)
      for (int x = 0; x < gx.block; ++x)
        ws.block()[std::size_t(y)*gx.block+x] = finite(finite(float(sources[0].row_ptr(reflect(oy+y,gy))[reflect(ox+x,gx)]))-base);
    fft.forward(ws.block().data(), ws.spectrum().data());
    fft.inverse(ws.spectrum().data(), ws.inverse().data());
    for (int y = 0; y < gy.block; ++y)
      for (int x = 0; x < gx.block; ++x)
        accum.row_ptr(oy+y)[ox+x] = ws.inverse()[std::size_t(y)*gx.block+x];
  } else if (algorithm == Algorithm::FFT3D) {
    const int c = T_slots / 2;
    NoisePower noise;
    if (denoise_) {
      if (parameters.primary_mode == PrimaryMode::Table) { noise.mode = PrimaryMode::Table; noise.table = parameters.primary; noise.multiplier = float(T_slots); }
      else noise.uniform = ((float(T_slots) * sigma_eff_) * sigma_eff_) / norm_;
    }
    const float lower = (beta_ - 1.0f) / beta_;
    const float* wx_a = wx_.analysis.data();
    const float* wx_s = wx_.synthesis.data();
    const std::size_t spatial_bins = fft.bins();

    for (int by = 0; by < gy.count; ++by) {
      const int oy = by * gy.step;
      if (!row.empty())
        std::memset(row.data(), 0,
                    static_cast<std::size_t>(row.stride_bytes()) * static_cast<std::size_t>(row.height()));

      for (int bx_start = 0; bx_start < gx.count; ++bx_start) {
        const int ox = bx_start * gx.step;

        if (T_slots == 1) {
          float* blk = ws.block(0).data();
          for (int y = 0; y < gy.block; ++y) {
            const float* src_row = ws.padded(0).row_ptr(oy + y) + ox;
            float* blk_row = blk + y * gx.block;
            const float wy = wy_.analysis[y];
            spatial_.gather_fft3d(src_row, wx_a, wy, blk_row, gx.block);
          }

          fft.forward(blk, ws.spectrum(0).data());
          auto* spec_b = ws.spectrum(0).data();
          const float scale = grid_.empty() ? 0 : (mean_scale_ * spec_b[0].real()) / grid_[0].real();
          if (denoise_) {
            kernel_(spec_b, grid_.empty() ? nullptr : grid_.data(), fft.bins(), scale, parameters);
            if (parameters.primary_mode == PrimaryMode::Table) enhance(spec_b);
          } else enhance(spec_b);
          fft.inverse(spec_b, ws.inverse(0).data());

          const float* inv_b = ws.inverse(0).data();
          for (int y = 0; y < gy.block; ++y) {
            float* r_row = row.row_ptr(y) + ox;
            const float* inv_row = inv_b + y * gx.block;
            spatial_.scatter_fft3d_block(inv_row, wx_s, r_row, gx.block);
          }
        } else {
          const std::complex<float>* spectra_ptrs[5];
          for (int j = 0; j < T_slots; ++j) {
            float* blk = ws.block(0).data();
            for (int y = 0; y < gy.block; ++y) {
              const float* src_row = ws.padded(j).row_ptr(oy + y) + ox;
              float* blk_row = blk + y * gx.block;
              const float wy = wy_.analysis[y];
              spatial_.gather_fft3d(src_row, wx_a, wy, blk_row, gx.block);
            }
            std::complex<float>* spec_j = ws.spectrum(0).data() + j * spatial_bins;
            fft.forward(blk, spec_j);
            spectra_ptrs[j] = spec_j;
          }

          std::complex<float>* out_spectrum = ws.spectrum(0).data() + T_slots * spatial_bins;
          temporal_kernel_(spectra_ptrs, T_slots, c, spatial_bins,
                           mean_scale_, grid_.empty() ? nullptr : grid_.data(),
                           noise, lower, out_spectrum);

          float* inv_buf = ws.inverse(0).data();
          enhance(out_spectrum);
          fft.inverse(out_spectrum, inv_buf);

          for (int y = 0; y < gy.block; ++y) {
            float* r_row = row.row_ptr(y) + ox;
            const float* inv_row = inv_buf + y * gx.block;
            spatial_.scatter_fft3d_block(inv_row, wx_s, r_row, gx.block);
          }
        }
      }

      if (!row.empty()) {
        for (int y = 0; y < gy.block; ++y) {
          const float* r_ptr = row.row_ptr(y);
          float* a_ptr = accum.row_ptr(oy + y);
          const float wy = wy_.synthesis[y];
          spatial_.scatter_fft3d_row(r_ptr, wy, a_ptr, gx.cover);
        }
      }
    }
  } else { // DFTTest
    if (T_slots == 1) {
      const float* h_data = h_.data();
      const float* h_syn_data = h_synthesis_.data();
      for (int by = 0; by < gy.count; ++by) {
        const int oy = by * gy.step;
        for (int bx_start = 0; bx_start < gx.count; ++bx_start) {
          const int ox = bx_start * gx.step;
          float* blk = ws.block(0).data();
          for (int y = 0; y < gy.block; ++y) {
            const float* src_row = ws.padded(0).row_ptr(oy + y) + ox;
            float* blk_row = blk + y * gx.block;
            const float* h_row = h_data + y * gx.block;
            spatial_.gather_dfttest(src_row, h_row, blk_row, gx.block);
          }

          fft.forward(blk, ws.spectrum(0).data());
          auto* spec_b = ws.spectrum(0).data();
          const float scale = grid_.empty() ? 0 : (mean_scale_ * spec_b[0].real()) / grid_[0].real();
          kernel_(spec_b, grid_.empty() ? nullptr : grid_.data(), fft.bins(), scale, parameters);
          fft.inverse(spec_b, ws.inverse(0).data());

          if (center_) {
            const int cy = gy.block / 2, cx = gx.block / 2;
            const auto ci = std::size_t(cy) * gx.block + cx;
            accum.row_ptr(oy + cy)[ox + cx] = ws.inverse(0).data()[ci] * h_syn_data[ci];
          } else {
            const float* inv_b = ws.inverse(0).data();
            for (int y = 0; y < gy.block; ++y) {
              float* acc_row = accum.row_ptr(oy + y) + ox;
              const float* inv_row = inv_b + y * gx.block;
              const float* h_syn_row = h_syn_data + y * gx.block;
              spatial_.scatter_dfttest(inv_row, h_syn_row, acc_row, gx.block);
            }
          }
        }
      }
    } else {
      const int c = T_slots / 2;
      const std::size_t bins_3d = std::size_t(T_slots) * gy.block * (gx.block / 2 + 1);
      const int spatial_block_size = gy.block * gx.block;

      for (int by = 0; by < gy.count; ++by) {
        const int oy = by * gy.step;
        for (int bx_start = 0; bx_start < gx.count; ++bx_start) {
          const int ox = bx_start * gx.step;

          float* block_3d = ws.block(0).data();
          for (int z = 0; z < T_slots; ++z) {
            const auto pad_z = ws.padded(z);
            const float* h_z = h_.data() + z * spatial_block_size;
            float* blk_z = block_3d + z * spatial_block_size;
            for (int y = 0; y < gy.block; ++y) {
              const float* src_row = pad_z.row_ptr(oy + y) + ox;
              float* blk_row = blk_z + y * gx.block;
              const float* h_row = h_z + y * gx.block;
              spatial_.gather_dfttest(src_row, h_row, blk_row, gx.block);
            }
          }

          std::complex<float>* spec_3d = ws.spectrum(0).data();
          fft3d_->forward(block_3d, spec_3d);

          const float g_ratio = (!grid_.empty() && mean_scale_ != 0.0f && grid_[0].real() != 0.0f)
                                    ? (mean_scale_ * spec_3d[0].real() / grid_[0].real())
                                    : 0.0f;
          kernel_(spec_3d, grid_.empty() ? nullptr : grid_.data(), bins_3d, g_ratio, parameters);

          float* inv_3d = ws.inverse(0).data();
          fft3d_->inverse(spec_3d, inv_3d);

          const float* inv_c = inv_3d + c * spatial_block_size;
          const float* h_syn_c = h_synthesis_.data() + c * spatial_block_size;

          if (center_) {
            const int cy = gy.block / 2, cx = gx.block / 2;
            const auto ci = std::size_t(cy) * gx.block + cx;
            accum.row_ptr(oy + cy)[ox + cx] = inv_c[ci] * h_syn_c[ci];
          } else {
            for (int y = 0; y < gy.block; ++y) {
              float* acc_row = accum.row_ptr(oy + y) + ox;
              const float* inv_row = inv_c + y * gx.block;
              const float* h_syn_row = h_syn_c + y * gx.block;
              spatial_.scatter_dfttest(inv_row, h_syn_row, acc_row, gx.block);
            }
          }
        }
      }
    }
  }

  if constexpr (std::is_same_v<T, float>) {
    const float scale = algorithm == Algorithm::FFT3D ? 1.0f : 1.0f / 255.0f;
    for (int y = 0; y < dst.height(); ++y) {
      const float* a_ptr = accum.row_ptr(y + gy.offset) + gx.offset;
      float* dst_row = dst.row_ptr(y);
      spatial_.store_output_float(a_ptr, dst_row, dst.width(), algorithm == Algorithm::FFT3D, scale);
    }
  } else if constexpr (sizeof(T) == 1) {
    const float peak = float((1 << format.bits) - 1);
    const float scale = float(1 << (format.bits - 8));
    for (int y = 0; y < dst.height(); ++y) {
      const float* a_ptr = accum.row_ptr(y + gy.offset) + gx.offset;
      auto* dst_row = reinterpret_cast<std::uint8_t*>(dst.row_ptr(y));
      spatial_.store_output_uint8(a_ptr, dst_row, dst.width(), algorithm == Algorithm::FFT3D, base, scale, peak);
    }
  } else {
    const float peak = float((1 << format.bits) - 1);
    const float scale = float(1 << (format.bits - 8));
    for (int y = 0; y < dst.height(); ++y) {
      const float* a_ptr = accum.row_ptr(y + gy.offset) + gx.offset;
      auto* dst_row = reinterpret_cast<std::uint16_t*>(dst.row_ptr(y));
      spatial_.store_output_uint16(a_ptr, dst_row, dst.width(), algorithm == Algorithm::FFT3D, base, scale, peak);
    }
  }
}
void Plan::process(span2d::Plane<const std::uint8_t> s, span2d::Plane<std::uint8_t> d) const {
  const span2d::Plane<const std::uint8_t> arr[1] = {s};
  process(span2d::Span<const span2d::Plane<const std::uint8_t>>(arr, 1), d);
}
void Plan::process(span2d::Plane<const std::uint16_t> s, span2d::Plane<std::uint16_t> d) const {
  const span2d::Plane<const std::uint16_t> arr[1] = {s};
  process(span2d::Span<const span2d::Plane<const std::uint16_t>>(arr, 1), d);
}
void Plan::process(span2d::Plane<const float> s, span2d::Plane<float> d) const {
  const span2d::Plane<const float> arr[1] = {s};
  process(span2d::Span<const span2d::Plane<const float>>(arr, 1), d);
}
void Plan::process(span2d::Plane<const std::uint8_t> s, span2d::Plane<std::uint8_t> d, runtime::Workspace& ws) const {
  const span2d::Plane<const std::uint8_t> arr[1] = {s};
  process(span2d::Span<const span2d::Plane<const std::uint8_t>>(arr, 1), d, ws);
}
void Plan::process(span2d::Plane<const std::uint16_t> s, span2d::Plane<std::uint16_t> d, runtime::Workspace& ws) const {
  const span2d::Plane<const std::uint16_t> arr[1] = {s};
  process(span2d::Span<const span2d::Plane<const std::uint16_t>>(arr, 1), d, ws);
}
void Plan::process(span2d::Plane<const float> s, span2d::Plane<float> d, runtime::Workspace& ws) const {
  const span2d::Plane<const float> arr[1] = {s};
  process(span2d::Span<const span2d::Plane<const float>>(arr, 1), d, ws);
}

void Plan::process(span2d::Span<const span2d::Plane<const std::uint8_t>> sources, span2d::Plane<std::uint8_t> dst) const {
  auto lease = pool_.acquire();
  run(sources, dst, *lease);
}
void Plan::process(span2d::Span<const span2d::Plane<const std::uint16_t>> sources, span2d::Plane<std::uint16_t> dst) const {
  auto lease = pool_.acquire();
  run(sources, dst, *lease);
}
void Plan::process(span2d::Span<const span2d::Plane<const float>> sources, span2d::Plane<float> dst) const {
  auto lease = pool_.acquire();
  run(sources, dst, *lease);
}

void Plan::process(span2d::Span<const span2d::Plane<const std::uint8_t>> sources, span2d::Plane<std::uint8_t> dst,
                   runtime::Workspace& ws) const {
  run(sources, dst, ws);
}
void Plan::process(span2d::Span<const span2d::Plane<const std::uint16_t>> sources, span2d::Plane<std::uint16_t> dst,
                   runtime::Workspace& ws) const {
  run(sources, dst, ws);
}
void Plan::process(span2d::Span<const span2d::Plane<const float>> sources, span2d::Plane<float> dst,
                   runtime::Workspace& ws) const {
  run(sources, dst, ws);
}
} // namespace neo_fft
