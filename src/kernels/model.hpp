#pragma once
#include <complex>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace neo_fft {
enum class SampleStorage { U8, U16, F32 };
template<class T> constexpr SampleStorage sample_storage = std::is_same_v<T,float> ? SampleStorage::F32 :
    sizeof(T) == 1 ? SampleStorage::U8 : SampleStorage::U16;

// Logical counts only, no padding. All operations check consumed intermediates.
// Buffers are disjoint except the explicitly in-place window/scale operations.
// On failure the private destination may be partially written, never published.
struct ModelKernels {
  void (*decode)(const void*, SampleStorage, float*, std::size_t, float base, float scale);
  void (*window)(float*, const float*, std::size_t, float first_factor);
  void (*power)(const std::complex<float>*, const std::complex<float>*, float ratio,
                float*, std::size_t, bool accumulate);
  // Ordered scalar addition is intentional: preserves first-minimum selection.
  float (*score)(const float*, const float*, std::size_t);
  void (*scale)(float*, const float*, std::size_t, float factor, float check_multiplier);
};
const ModelKernels& model_scalar();
ModelKernels select_model(int opt);
} // namespace neo_fft
