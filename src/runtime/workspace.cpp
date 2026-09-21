#include "runtime/workspace.hpp"
#include "base/checked.hpp"
#include <cstring>

namespace neo_fft::runtime {

WorkspaceBudget make_workspace_budget(const Geometry& geom, std::size_t samples, std::size_t bins, bool has_row_buffer,
                                      int temporal_slots, int batch_size) {
  require(batch_size > 0, "batch_size must be positive");
  require(temporal_slots > 0, "temporal_slots must be positive");
  WorkspaceBudget b;
  b.batch_size = batch_size;
  b.temporal_slots = temporal_slots;
  b.samples = samples;
  b.bins = bins;

  b.accum_width = geom.x.cover;
  b.accum_height = geom.y.cover;
  b.accum_stride_bytes = static_cast<std::ptrdiff_t>(
      align_up(mul_size(static_cast<std::size_t>(b.accum_width), sizeof(float)), kSimdAlignment));
  const auto accum_bytes =
      align_up(mul_size(static_cast<std::size_t>(b.accum_stride_bytes), static_cast<std::size_t>(b.accum_height)),
               kSimdAlignment);

  b.padded_width = geom.x.cover;
  b.padded_height = geom.y.cover;
  b.padded_stride_bytes = static_cast<std::ptrdiff_t>(
      align_up(mul_size(static_cast<std::size_t>(b.padded_width), sizeof(float)), kSimdAlignment));
  const auto single_padded_bytes =
      mul_size(static_cast<std::size_t>(b.padded_stride_bytes), static_cast<std::size_t>(b.padded_height));
  b.padded_bytes = align_up(mul_size(single_padded_bytes, static_cast<std::size_t>(temporal_slots)), kSimdAlignment);

  std::size_t row_bytes = 0;
  if (has_row_buffer) {
    b.row_width = geom.x.cover;
    b.row_height = geom.y.block;
    b.row_stride_bytes = static_cast<std::ptrdiff_t>(
        align_up(mul_size(static_cast<std::size_t>(b.row_width), sizeof(float)), kSimdAlignment));
    row_bytes =
        align_up(mul_size(static_cast<std::size_t>(b.row_stride_bytes), static_cast<std::size_t>(b.row_height)),
                 kSimdAlignment);
  }

  const auto block_bytes =
      align_up(mul_size(mul_size(static_cast<std::size_t>(batch_size), b.samples), sizeof(float)), kSimdAlignment);
  const auto inverse_bytes =
      align_up(mul_size(mul_size(static_cast<std::size_t>(batch_size), b.samples), sizeof(float)), kSimdAlignment);
  const auto spectrum_bytes = align_up(
      mul_size(mul_size(static_cast<std::size_t>(batch_size), b.bins), sizeof(std::complex<float>)), kSimdAlignment);

  b.accum_offset = 0;
  b.padded_offset = b.accum_offset + accum_bytes;
  b.row_offset = b.padded_offset + b.padded_bytes;
  b.block_offset = b.row_offset + row_bytes;
  b.inverse_offset = b.block_offset + block_bytes;
  b.spectrum_offset = b.inverse_offset + inverse_bytes;
  b.total_bytes = b.spectrum_offset + spectrum_bytes;

  return b;
}

WorkspaceBudget make_workspace_budget(const Geometry& geom, const RealFFT& fft, bool has_row_buffer, int batch_size) {
  return make_workspace_budget(geom, fft.samples(), fft.bins(), has_row_buffer, 1, batch_size);
}

Workspace::Workspace(const WorkspaceBudget& budget)
    : budget_(budget), storage_(budget.total_bytes + kSimdAlignment) {
  const auto raw = reinterpret_cast<std::uintptr_t>(storage_.data());
  const auto aligned = (raw + (kSimdAlignment - 1)) & ~(std::uintptr_t(kSimdAlignment - 1));
  auto* base = reinterpret_cast<std::byte*>(aligned);

  accum_ptr_ = reinterpret_cast<float*>(base + budget_.accum_offset);
  padded_ptr_ = reinterpret_cast<float*>(base + budget_.padded_offset);
  if (budget_.row_height > 0) {
    row_ptr_ = reinterpret_cast<float*>(base + budget_.row_offset);
  }
  block_ptr_ = reinterpret_cast<float*>(base + budget_.block_offset);
  inverse_ptr_ = reinterpret_cast<float*>(base + budget_.inverse_offset);
  spectrum_ptr_ = reinterpret_cast<std::complex<float>*>(base + budget_.spectrum_offset);
}

void Workspace::reset() noexcept {
  if (accum_ptr_ && budget_.accum_height > 0 && budget_.accum_stride_bytes > 0) {
    std::memset(accum_ptr_, 0,
                static_cast<std::size_t>(budget_.accum_stride_bytes) * static_cast<std::size_t>(budget_.accum_height));
  }
}

} // namespace neo_fft::runtime
