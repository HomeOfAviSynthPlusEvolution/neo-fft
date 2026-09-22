#pragma once
#include <dualsynth/span2d.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace neo_fft {
// Literal diagnostics must not allocate on the successful hot path.
inline void require(bool ok, const char* message) {
  if (!ok)
    throw std::invalid_argument(message);
}
inline void require(bool ok, const std::string& message) {
  if (!ok)
    throw std::invalid_argument(message);
}
inline std::size_t add_size(std::size_t a, std::size_t b) {
  require(b <= std::numeric_limits<std::size_t>::max() - a, "storage addition overflow");
  return a + b;
}
inline std::size_t mul_size(std::size_t a, std::size_t b) {
  require(a == 0 || b <= std::numeric_limits<std::size_t>::max() / a, "storage product overflow");
  return a * b;
}
inline int dimension(std::int64_t n) {
  require(n > 0 && n <= INT32_MAX, "dimension outside int32 domain");
  return static_cast<int>(n);
}
inline std::int64_t ceil_div(std::int64_t n, std::int64_t d) {
  require(n >= 0 && d > 0, "invalid ceiling division");
  return n / d + (n % d != 0);
}
inline float finite(float v) {
  if (!std::isfinite(v))
    throw std::runtime_error("non-finite sample or intermediate");
  return v;
}
template <class T>
std::size_t plane_extent(std::int64_t w, std::int64_t h, std::ptrdiff_t stride) {
  dimension(w);
  dimension(h);
  require(stride > 0 && stride % sizeof(T) == 0 && stride / sizeof(T) <= INT32_MAX, "invalid plane byte stride");
  const auto row = mul_size(static_cast<std::size_t>(w), sizeof(T));
  require(static_cast<std::size_t>(stride) >= row, "stride shorter than active row");
  const auto extent = add_size(mul_size(static_cast<std::size_t>(h - 1), static_cast<std::size_t>(stride)), row);
  require(extent <= static_cast<std::size_t>(PTRDIFF_MAX), "plane extent exceeds ptrdiff_t");
  return extent;
}
template <class T>
span2d::Plane<T> checked_plane(T* data, std::int64_t w, std::int64_t h, std::ptrdiff_t stride, std::size_t accessible) {
  const auto extent = plane_extent<T>(w, h, stride);
  require(data && reinterpret_cast<std::uintptr_t>(data) % alignof(T) == 0, "invalid plane pointer alignment");
  require(extent <= accessible, "plane exceeds owner extent");
  require(extent <= UINTPTR_MAX - reinterpret_cast<std::uintptr_t>(data), "plane address overflow");
  return {data, static_cast<int>(w), static_cast<int>(h), stride};
}
template <class T>
span2d::Plane<T> checked_subplane(span2d::Plane<T> p, int x, int y, int w, int h) {
  require(x >= 0 && y >= 0 && w > 0 && h > 0 && std::int64_t(x) + w <= p.width() && std::int64_t(y) + h <= p.height(),
          "subplane outside active extent");
  return p.subplane(x, y, w, h);
}
template <class T>
std::vector<T> buffer(std::size_t count) {
  const auto bytes = mul_size(count, sizeof(T));
  require(bytes <= static_cast<std::size_t>(PTRDIFF_MAX) && count <= std::vector<T>().max_size(),
          "allocation extent unrepresentable");
  return std::vector<T>(count); // Constructs live T objects, including complex storage.
}
template <class T>
void disjoint(const T* a, std::size_t an, const void* b, std::size_t bn) {
  const auto ap = reinterpret_cast<std::uintptr_t>(a), bp = reinterpret_cast<std::uintptr_t>(b);
  require(an <= UINTPTR_MAX - ap && bn <= UINTPTR_MAX - bp, "address extent overflow");
  require(ap + an <= bp || bp + bn <= ap, "partially overlapping or aliased buffers");
}
} // namespace neo_fft
