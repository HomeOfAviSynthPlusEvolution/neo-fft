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
  require(c.bt >= 1 && c.bt <= 5, "FFT3D bt outside 1..5");
  nonnegative(c.sigma, "FFT3D sigma");
  nonnegative(c.degrid, "FFT3D degrid");
  require(std::isfinite(c.beta) && c.beta >= 1, "FFT3D beta must be >=1");
  select_spectral(c.opt);
}
void validate(const DFTConfig& c) {
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
    : geometry(geometry3d(w, h, c)), format(f), algorithm(Algorithm::FFT3D), temporal_size(c.bt), fft(c.bh, c.bw),
      wx_(fft3d_window(c.bw, geometry.x.overlap, c.wintype)), wy_(fft3d_window(c.bh, geometry.y.overlap, c.wintype)),
      kernel_(select_spectral(c.opt)), spatial_(select_spatial(c.opt)), mean_scale_(c.degrid),
      pool_(runtime::make_workspace_budget(geometry, fft.samples(), fft.bins() * std::size_t(std::max(1, c.bt)), true,
                                           std::max(1, c.bt),
                                           select_optimal_batch_size(fft.samples(), geometry.x.count))) {
  valid_format(f);
  const float factor = f.floating ? 1.0f / 255 : float(1 << (f.bits - 8));
  sigma_eff_ = finite(c.sigma * factor);
  norm_ = 1.0f / (float(c.bw) * float(c.bh));
  beta_ = c.beta;
  params_.a = finite((sigma_eff_ * sigma_eff_) / norm_);
  params_.floor = (c.beta - 1) / c.beta;
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
Plan::Plan(int w, int h, SampleFormat f, const DFTConfig& c)
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
  const float scale = c.ftype < 2 ? win.wscale : 1.0f;
  params_ = {c.ftype,
             finite(c.sigma / scale),
             finite(c.sigma2 / scale),
             finite(c.pmin / win.wscale),
             finite(c.pmax / win.wscale),
             c.f0beta,
             0};
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
template <class T>
void Plan::run(span2d::Span<const span2d::Plane<const T>> sources, span2d::Plane<T> dst, runtime::Workspace& ws) const {
  require(!sources.empty(), "sources must not be empty");
  const auto& gx = geometry.x;
  const auto& gy = geometry.y;
  const int T_slots = int(sources.size());
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
      for (int y = 0; y < src.height(); ++y)
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

  ws.reset();
  for (int j = 0; j < T_slots; ++j) {
    pad_source(sources[j], ws.padded(j), geometry, format, algorithm);
  }

  auto accum = ws.accum();
  auto row = algorithm == Algorithm::FFT3D ? ws.row() : span2d::Plane<float>{};
  const float base =
      !format.floating && format.chroma && algorithm == Algorithm::FFT3D ? float(1 << (format.bits - 1)) : 0;

  if (algorithm == Algorithm::FFT3D) {
    const int c = T_slots / 2;
    const float noise = ((float(T_slots) * sigma_eff_) * sigma_eff_) / norm_;
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
          kernel_(spec_b, grid_.empty() ? nullptr : grid_.data(), fft.bins(), scale, params_);
          fft.inverse(spec_b, ws.inverse(0).data());

          const float* inv_b = ws.inverse(0).data();
          for (int y = 0; y < gy.block; ++y) {
            float* r_row = row.row_ptr(y) + ox;
            const float* inv_row = inv_b + y * gx.block;
            spatial_.scatter_fft3d_block(inv_row, wx_s, r_row, gx.block);
          }
        } else {
          std::vector<std::vector<std::complex<float>>> spectra_slots(T_slots, std::vector<std::complex<float>>(spatial_bins));
          std::vector<const std::complex<float>*> spectra_ptrs(T_slots);
          std::vector<float> block_buf(fft.samples());

          for (int j = 0; j < T_slots; ++j) {
            for (int y = 0; y < gy.block; ++y) {
              const float* src_row = ws.padded(j).row_ptr(oy + y) + ox;
              float* blk_row = block_buf.data() + y * gx.block;
              const float wy = wy_.analysis[y];
              spatial_.gather_fft3d(src_row, wx_a, wy, blk_row, gx.block);
            }
            fft.forward(block_buf.data(), spectra_slots[j].data());
            spectra_ptrs[j] = spectra_slots[j].data();
          }

          std::vector<std::complex<float>> out_spectrum(spatial_bins);
          fft3d_temporal_filter(spectra_ptrs.data(), T_slots, c, spatial_bins,
                               mean_scale_, grid_.empty() ? nullptr : grid_.data(),
                               noise, lower, out_spectrum.data());

          std::vector<float> inv_buf(fft.samples());
          fft.inverse(out_spectrum.data(), inv_buf.data());

          for (int y = 0; y < gy.block; ++y) {
            float* r_row = row.row_ptr(y) + ox;
            const float* inv_row = inv_buf.data() + y * gx.block;
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
          kernel_(spec_b, grid_.empty() ? nullptr : grid_.data(), fft.bins(), scale, params_);
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
      const std::size_t volume_3d = std::size_t(T_slots) * gy.block * gx.block;
      const std::size_t bins_3d = std::size_t(T_slots) * gy.block * (gx.block / 2 + 1);
      const int spatial_block_size = gy.block * gx.block;

      for (int by = 0; by < gy.count; ++by) {
        const int oy = by * gy.step;
        for (int bx_start = 0; bx_start < gx.count; ++bx_start) {
          const int ox = bx_start * gx.step;

          std::vector<float> block_3d(volume_3d);
          for (int z = 0; z < T_slots; ++z) {
            const auto pad_z = ws.padded(z);
            const float* h_z = h_.data() + z * spatial_block_size;
            float* blk_z = block_3d.data() + z * spatial_block_size;
            for (int y = 0; y < gy.block; ++y) {
              const float* src_row = pad_z.row_ptr(oy + y) + ox;
              float* blk_row = blk_z + y * gx.block;
              const float* h_row = h_z + y * gx.block;
              spatial_.gather_dfttest(src_row, h_row, blk_row, gx.block);
            }
          }

          std::vector<std::complex<float>> spec_3d(bins_3d);
          fft3d_->forward(block_3d.data(), spec_3d.data());

          const float g_ratio = (!grid_.empty() && mean_scale_ != 0.0f && grid_[0].real() != 0.0f)
                                    ? (mean_scale_ * spec_3d[0].real() / grid_[0].real())
                                    : 0.0f;
          kernel_(spec_3d.data(), grid_.empty() ? nullptr : grid_.data(), bins_3d, g_ratio, params_);

          std::vector<float> inv_3d(volume_3d);
          fft3d_->inverse(spec_3d.data(), inv_3d.data());

          const float* inv_c = inv_3d.data() + c * spatial_block_size;
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
