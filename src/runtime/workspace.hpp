#pragma once
#include "spectral/geometry.hpp"
#include "spectral/fft.hpp"
#include <dualsynth/span2d.hpp>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace neo_fft::runtime {

inline constexpr std::size_t kSimdAlignment = 64;

inline std::size_t align_up(std::size_t n, std::size_t alignment = kSimdAlignment) {
  require(alignment && !(alignment & (alignment-1)), "alignment must be a power of two");
  return add_size(n, alignment - 1) & ~(alignment - 1);
}

struct WorkspaceBudget {
  int accum_width = 0;
  int accum_height = 0;
  std::ptrdiff_t accum_stride_bytes = 0;

  int padded_width = 0;
  int padded_height = 0;
  std::ptrdiff_t padded_stride_bytes = 0;

  int row_width = 0;
  int row_height = 0;
  std::ptrdiff_t row_stride_bytes = 0;

  std::size_t samples = 0;
  std::size_t bins = 0;
  int batch_size = 1;
  int temporal_slots = 1;

  std::size_t accum_offset = 0;
  std::size_t padded_offset = 0;
  std::size_t row_offset = 0;
  std::size_t block_offset = 0;
  std::size_t inverse_offset = 0;
  std::size_t spectrum_offset = 0;
  std::size_t padded_bytes = 0;
  std::size_t total_bytes = 0;
};

WorkspaceBudget make_workspace_budget(const Geometry& geom, std::size_t samples, std::size_t bins, bool has_row_buffer,
                                      int temporal_slots = 1, int batch_size = 1);
WorkspaceBudget make_workspace_budget(const Geometry& geom, const RealFFT& fft, bool has_row_buffer, int batch_size = 1);

class Workspace {
public:
  explicit Workspace(const WorkspaceBudget& budget);
  ~Workspace() = default;

  Workspace(const Workspace&) = delete;
  Workspace& operator=(const Workspace&) = delete;
  Workspace(Workspace&&) noexcept = default;
  Workspace& operator=(Workspace&&) noexcept = default;

  const WorkspaceBudget& budget() const noexcept { return budget_; }

  // 2D planes
  span2d::Plane<float> accum() const noexcept {
    return span2d::Plane<float>(accum_ptr_, budget_.accum_width, budget_.accum_height, budget_.accum_stride_bytes);
  }

  span2d::Plane<float> padded(int slot = 0) const noexcept {
    const std::size_t plane_floats = static_cast<std::size_t>(budget_.padded_height) *
                                     (static_cast<std::size_t>(budget_.padded_stride_bytes) / sizeof(float));
    return span2d::Plane<float>(padded_ptr_ + slot * plane_floats, budget_.padded_width, budget_.padded_height,
                                budget_.padded_stride_bytes);
  }

  span2d::Plane<float> row() const noexcept {
    return span2d::Plane<float>(row_ptr_, budget_.row_width, budget_.row_height, budget_.row_stride_bytes);
  }

  // 1D continuous spans
  span2d::Span<float> block(int batch_index = 0) const noexcept {
    return span2d::Span<float>(block_ptr_ + batch_index * budget_.samples, budget_.samples);
  }

  span2d::Span<float> inverse(int batch_index = 0) const noexcept {
    return span2d::Span<float>(inverse_ptr_ + batch_index * budget_.samples, budget_.samples);
  }

  span2d::Span<std::complex<float>> spectrum(int batch_index = 0) const noexcept {
    return span2d::Span<std::complex<float>>(spectrum_ptr_ + batch_index * budget_.bins, budget_.bins);
  }

  // Restrict views for compiler vectorization without aliasing overhead
  span2d::RestrictPlane<float> accum_restrict() const noexcept { return accum().as_restrict(); }
  span2d::RestrictPlane<float> padded_restrict() const noexcept { return padded().as_restrict(); }
  span2d::RestrictPlane<float> row_restrict() const noexcept { return row().as_restrict(); }
  span2d::RestrictSpan<float> block_restrict(int batch_index = 0) const noexcept {
    return block(batch_index).as_restrict();
  }
  span2d::RestrictSpan<float> inverse_restrict(int batch_index = 0) const noexcept {
    return inverse(batch_index).as_restrict();
  }

  // Zero accumulator plane before processing a new frame/plane
  std::size_t retained_bytes() const { return add_size(sizeof(*this),storage_.capacity()); }
  void reset() noexcept;

private:
  WorkspaceBudget budget_;
  std::vector<std::byte> storage_;
  float* accum_ptr_ = nullptr;
  float* padded_ptr_ = nullptr;
  float* row_ptr_ = nullptr;
  float* block_ptr_ = nullptr;
  float* inverse_ptr_ = nullptr;
  std::complex<float>* spectrum_ptr_ = nullptr;
};

} // namespace neo_fft::runtime
