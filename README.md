# neo-fft

**English** | [简体中文](README.zh-CN.md) | [日本語](README.ja.md)

neo-fft is a frequency-domain filtering plugin for VapourSynth and AviSynth. It independently rewrites Neo FFT3D and Neo DFTTest as a single plugin, providing spatial and temporal denoising, Kalman filtering, frequency-domain sharpening, dehalo, and configurable spectral filtering.

The implementation uses C++17, with scalar kernels, cross-platform SIMD through Google Highway, and PocketFFT-based transforms. DualSynth2 connects the computation core to both hosts. VapourSynth uses `core.neo_fft`; AviSynth uses functions prefixed with `neo_fft_`.

## Design

neo-fft separates frequency-domain calculations from host frame management. The core handles block extraction, windows, FFT, spectral filtering, and reconstruction; the host layer manages parameters, frame requests, properties, and output allocation. FFT3D and DFTTest share infrastructure while retaining their own window, noise-model, and overlap-add rules. The core can be built and tested independently.

The implementation was developed from behavioral specifications. Scalar implementations and independent mathematical definitions are used to check calculations, followed by verification of SIMD paths and public behavior against fixed reference plugins. Identical parameters do not guarantee pixel-identical output to every historical FFT3DFilter, DFTTest, or Neo build. Floating-point rounding, threshold branches, Kalman state handling, and dithering can affect results. See the [migration guide](docs/api/en/migration.md) for specific differences.

Both filters and FFT execute on the calling host thread, without creating worker threads or thread pools. The host can still request multiple frames concurrently; SIMD and batched FFT do not imply internal multithreading.

## Supported operations

| Function | Purpose |
|---|---|
| `FFT3D` | Single-frame or 2–5-frame Wiener denoising, Kalman filtering, sharpening, and dehalo, with frequency-dependent and sampled noise models. |
| `DFTTest` | Five spectral filter types, spatial and temporal overlap-add, twelve windows, frequency curves, noise sampling, and dithering for 8-bit output. |
| `KernelInfo` | Automatically selected SIMD target, FFT implementation, and vector width. |

Both video filters support constant-format, constant-size planar GRAY/YUV/RGB with 8/10/12/14/16-bit integer or 32-bit floating-point samples. AviSynth also supports planar YUVA/RGBA. Output dimensions, format, frame count, and frame rate match the input. Block geometry and boundary requirements are checked for each processed plane.

All non-alpha planes are processed by default; alpha is copied. `planes=[0]` processes only the first plane, while `planes=[]` processes none. AviSynth also accepts the legacy `y/u/v/a` plane modes. Mode 1 explicitly leaves the output plane unwritten, for scripts that subsequently discard that plane. See the API reference for precedence rules.

FFT3D defaults to `bt=3`, using three neighboring frames; `bt=1` selects spatial-only denoising, and `bt=0` selects Kalman filtering. Kalman uses at most eight historical warmup frames by default and can reuse nearby checkpoints, avoiding replay from the beginning of the clip after a seek. Cache contents and request history can affect the recursive result. The raw-spectrum cache for ordinary temporal filtering defaults to a 128 MiB budget, adjustable through `cache_mb` and `cache_frames`; this is not a process-wide memory limit.

DFTTest defaults to `tbsize=1`, using only the current frame. Increasing `tbsize` enables temporal filtering, and `tmode` selects center output or temporal overlap-add. Neither filter estimates motion vectors or performs motion compensation.

## Documentation and use

The API reference explains how to call the functions. The knowledge base explains how windows, transforms, spectral models, and reconstruction turn inputs into outputs.

- [API reference](docs/api/en/README.md): signatures, parameters, defaults, and usage examples.
- [Knowledge base](docs/knowledge/en/README.md): representations, formulas, operation order, boundaries, and precision.
- [Migrating from legacy Neo filters](docs/api/en/migration.md): required script changes and output differences.

Load the built plugin explicitly, or place it in VapourSynth's plugin autoload directory. The example below uses the Windows filename; Linux uses `neo-fft.so`. Substitute the actual plugin path for your platform. The VapourSynth plugin identifier is `org.neofilters.neo_fft`.

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")

print(core.neo_fft.KernelInfo())
clip = core.std.BlankClip(width=640, height=360, format=vs.YUV420P8, length=24)
output = core.neo_fft.FFT3D(clip, sigma=2.0, bt=3, planes=[0])
# Alternative: apply DFTTest spatial denoising to the original input.
# output = core.neo_fft.DFTTest(clip, sigma=8.0, tbsize=1, planes=[0])
output.set_output()
```

This minimal example uses a synthetic clip to demonstrate the calls. `sigma` has different mathematical meanings in the two filters: matching its value does not match denoising strength. Adjust each filter separately.

The same plugin also provides an AviSynth C++ interface. Use `LoadPlugin` with a host that supports interface version 11:

```avs
LoadPlugin("/path/to/neo-fft.dll")
clip = BlankClip(width=640, height=360, length=24, pixel_type="YV12")
return neo_fft_FFT3D(clip, sigma=2.0, bt=3, planes=[0])
# Alternative DFTTest call:
# return neo_fft_DFTTest(clip, sigma=8.0, tbsize=1, planes=[0])
```

AviSynth array parameters accept native arrays such as `[0, 1]`; a scalar is shorthand for one element. DFTTest curves and noise-sampling locations also accept numeric strings separated by whitespace, commas, or colons. Both filters forward input audio and parity, and output frame properties come from the corresponding source frame.

To explicitly express an omitted parameter, use `None` in VapourSynth or `Undefined()` in AviSynth. For example, `planes=None` or `planes=Undefined()` uses the default plane selection; the latter also lets AVS `y/u/v/a` take effect. An empty array `[]` explicitly selects no planes.

FFT3D's `mt/ncpu/measure`, DFTTest's `threads/fft_threads`, and both filters' `fft_backend` are accepted only for compatibility and completely ignored. They neither create threads nor select an FFT backend. Named arguments are recommended for migrating old scripts; see the API reference for parameter order and migration restrictions.

## SIMD and CPU selection

SIMD builds select a compiled Highway target supported by the running CPU. Scalar fallback remains available. FFT has independent implementation dispatch, with specialized transform paths for some fixed sizes.

For both filters, `opt=1` selects scalar kernels owned by this project, but does not force scalar FFT. Other supported values select automatic SIMD; legacy `opt` values cannot select or limit an ISA. Use `NEO_FFT_ENABLE_SIMD=OFF` to disable SIMD throughout the build.

`core.neo_fft.KernelInfo()` returns a dictionary containing `fft_backend`, `target`, `fft`, and `fft_lanes`, accessed by field name. AviSynth's `neo_fft_KernelInfo()` returns an array in the fixed order `[fft_backend, target, fft, fft_lanes]`. The query describes automatic selection, not an execution trace of a specific filter instance. `fft_lanes` is a vector lane count, not a thread count or a speedup estimate.

The FFT target can differ from the general-kernel target. Wider SIMD does not guarantee higher throughput. See [KernelInfo](docs/knowledge/en/kernel-info.md) and [execution and precision](docs/knowledge/en/shared/execution-precision.md).

## Building and testing

Requires CMake 3.24 or later, Git, and a C++17 compiler. CMake retrieves pinned DualSynth2 and PocketFFT sources, plus Highway 1.4.0 when SIMD is enabled. Both host SDKs are discovered locally or downloaded automatically.

The following builds the dual-host plugin and core tests without requiring an installed video host. On Windows, the default AVS C++ interface requires MSVC or clang-cl; MinGW builds must disable the AVS interface.

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DNEO_FFT_TEST_VAPOURSYNTH=OFF
cmake --build build/release --config Release --parallel 4
ctest --test-dir build/release -C Release --output-on-failure
```

| Option | Purpose |
|---|---|
| `NEO_FFT_BUILD_VAPOURSYNTH=OFF` | Disable the VapourSynth entry point. |
| `NEO_FFT_BUILD_AVISYNTH=OFF` | Disable the AviSynth entry point; setting both host options to OFF builds only the core. |
| `NEO_FFT_ENABLE_SIMD=OFF` | Disable SIMD kernels and vectorized FFT configurations. |
| `BUILD_TESTING=OFF` | Do not build tests. |
| `NEO_FFT_TEST_VAPOURSYNTH=OFF` | Keep core tests and the plugin, but skip VapourSynth host tests; defaults to ON when the VS interface and testing are enabled. |
| `Python3_EXECUTABLE=/path/to/python` | Select Python for host tests; VS tests require an environment that can import VapourSynth and load a matching-architecture runtime. |
| `NEO_FFT_VS_SDK=/path/to/sdk` | Specify a local VapourSynth SDK. |
| `NEO_FFT_AVS_SDK=/path/to/sdk` | Specify a local AviSynth SDK. |
| `NEO_FFT_TEST_AVISYNTH=ON` | Enable AviSynth host tests, disabled by default; requires Python and a matching-architecture runtime. |
| `NEO_FFT_AVISYNTH_RUNTIME=/path/to/avisynth.dll` | Specify the runtime library for AviSynth host tests. |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | Use local DualSynth2 sources instead of downloading the pinned version. |

The plugin target is `neo_fft`, and the output basename is `neo-fft`, with both host entry points enabled by default. Tests cover FFT, windows, filter operators, boundaries, scalar/SIMD comparisons, and host behavior. Black-box comparisons with old plugins additionally require fixed reference binaries.

CI covers Windows x64, Linux x64, macOS ARM64, and Linux ASan/UBSan checks. Linux runners use Ubuntu 26.04 with the system-default GCC; sanitizer checks use Clang 22. The release workflow builds Windows, Linux, and macOS x64/ARM64 artifacts. VapourSynth host tests currently run on Windows x64; AviSynth host tests are enabled separately through the options above. Packages record their actual build environment and test coverage; they are not necessarily compatible with every Linux distribution.

## Performance

Existing measurements put throughput on common FFT3D paths at approximately **1.72–3.05× that of legacy Neo FFT3D**, and common DFTTest Wiener paths at **4.58–7.04× that of legacy Neo DFTTest**. The ratio is **neo-fft throughput / reference-filter throughput**, equivalently reference-filter time / neo-fft time; **greater than 1 means neo-fft is faster**. The ranges below span 8/16-bit integer, 32-bit float, and AVX2/AVX-512 configurations; they are not confidence intervals.

| Common path | Relative throughput |
|---|---:|
| FFT3D spatial denoising (`bt=1`) | 2.39–2.85× |
| FFT3D two-frame denoising (`bt=2`) | 2.32–3.05× |
| FFT3D three-frame denoising (`bt=3`, default temporal mode) | 1.97–2.53× |
| FFT3D four-frame denoising (`bt=4`) | 1.86–2.37× |
| FFT3D five-frame denoising (`bt=5`) | 1.72–2.15× |
| FFT3D Kalman (`bt=0`, sequential requests) | 1.90–2.43× |
| DFTTest spatial Wiener (`tbsize=1`) | 4.60–6.65× |
| DFTTest three-frame Wiener (`tbsize=3, tmode=0`) | 5.25–7.04× |
| DFTTest five-frame Wiener (`tbsize=5, tmode=0`) | 4.58–5.77× |

These ranges come from historical single-threaded AviSynth comparisons against Highway-modernized versions of legacy Neo FFT3D / Neo DFTTest. The measurements exclude subsequent optimizations and do not represent the throughput of an entire processing chain. Actual results vary with input, parameters, hardware, and host concurrency.

## Development and contributions

Maintainers are responsible for technical direction, change review, and releases. Bug reports, suggestions, and contributions are welcome. Please discuss goals and approaches before changing numerical semantics, public interfaces, or major architecture.

This project uses AI assistance for implementation, testing, and review. Contributions should explain the problem, approach, validation, and how AI was involved. Reports should include the version, OS, CPU, compiler, build options, input/output formats, and a minimal reproducer. Numerical-difference reports should also identify the reference version, parameters, and request order; performance reports should describe the timing scope and thread configuration.

## Acknowledgments and license

Thanks to the authors and contributors of the following upstream projects, whose work provided the foundation for neo-fft's interface and frequency-domain filtering features:

- [FFT3DFilter](https://github.com/pinterf/fft3dfilter): originally developed by Alexander G. Balakhnin (Fizick), adapted for AviSynth 2.6 by martin53, and further improved with high-bit-depth support by Ferenc Pintér (pinterf).
- [DFTTest](https://github.com/pinterf/dfttest): originally developed by tritical, with 16-bit processing by Firesledge, an AviSynth+ port by DJATOM, and further high-bit-depth and cross-platform work by pinterf.

neo-fft also uses the following libraries:

- [Google Highway](https://github.com/google/highway): cross-platform SIMD support.
- [PocketFFT](https://github.com/mreineck/pocketfft): frequency-domain transforms for FFT3D and DFTTest.
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2): connects VapourSynth and AviSynth to the shared computation core.

Thanks to the developers and users who contribute tests, reports, and improvements.

Thanks to [SB.SB](https://sb.sb) for sponsoring the LLM subscription used in this project's development.

neo-fft is licensed under the GNU General Public License, version 2 or later (`GPL-2.0-or-later`). See [LICENSE](LICENSE). Third-party components retain their own copyright notices and license terms.
