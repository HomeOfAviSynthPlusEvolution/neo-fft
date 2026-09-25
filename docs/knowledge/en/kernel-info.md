# Interpreting kernel and FFT dispatch

[Contents](README.md)

KernelInfo describes automatic process-level selection, not a trace of a filter instance. See [API invocation and field types](../../api/en/kernel-info.md).

## Version and binary identity

The project VERSION in CMakeLists.txt is the single version source. Configuration generates the shared host header and Windows DLL version resource. Numeric VS registration stores only major and minor; the plugin description or DLL product version provides all three components. A version number does not distinguish local modifications under the same version. For reproducibility, also record the source commit, build configuration and binary checksum. KernelInfo SIMD/FFT fields cannot replace this identity information.

## Two independent selections

Own kernels use Highway runtime dispatch based on CPU/build capabilities. FFT uses a separate PocketFFT implementation configuration. opt=1 in either filter selects only own scalar operators. Other accepted opt values select automatic dispatch rather than a numerically named fixed ISA.

Thus target and fft_backend need not have matching names, and target cannot determine FFT lane count. `fft` identifies the selected profile, `fft_backend` the actual implementation, and `fft_lanes` its float SIMD width.

## Reading an example

Suppose a query reports `target=AVX3_SPR`, `fft=pocketfft-native`, `fft_backend=pocketfft-avx512`, `fft_lanes=16`. Automatic own kernels chose that Highway target; the native FFT profile selected AVX512 with 16 float lanes per vector. This does not mean 16 threads.

Creating an `opt=1` DFTTest in that process does not change the query's automatic target. The instance uses scalar own kernels while FFT may remain vectorized. Unsupported or unbuilt CPU targets fall back to available implementations rather than executing unsupported instructions.

## Specialized block paths

Some linear-Wiener DFTTest paths batch fixed FFT groups and some volume sizes have specialized transforms. Eligible tmode=0 paths can compute only the required inverse temporal center. Other types, exponents, sizes or OLA paths use their corresponding fallback. KernelInfo cannot prove that a call used a particular codelet or translate lane count directly into speedup.

Both filters and FFT are internally single-threaded. Compatibility inputs such as threads, fft_threads and fft_backend are completely ignored; KernelInfo reports no thread counts. The read-only fft_backend diagnostic reports the actual implementation independently of the ignored input. See [execution precision](shared/execution-precision.md).
