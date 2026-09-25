# neo-fft API reference

[Contents](README.md)

These pages cover public calling conventions, formats, types, defaults, return values, runnable examples and common errors. Formulas, processing order and numerical examples live in the [knowledge base](../../knowledge/en/README.md).

VapourSynth uses core.neo_fft, with plugin identifier org.neofilters.neo_fft. AviSynth uses the neo_fft_ function prefix. Both hosts expose FFT3D, DFTTest and KernelInfo.

| Function | Purpose and return value |
|---|---|
| [FFT3D](fft3d.md) | Spatial/temporal denoising, Kalman, sharpening and dehalo; VideoNode |
| [DFTTest](dfttest.md) | Five spectral filters, spatial/temporal overlap, curves/sampling and dither; VideoNode |
| [KernelInfo](kernel-info.md) | Automatic SIMD and FFT configuration; diagnostic dictionary |

For older Neo FFT3D / Neo DFTTest scripts, start with the [migration guide](migration.md), which separates required edits from differences in results.

## Version identification

The first release version is `0.9.0`. The VS plugin description and AVS loading description contain the full version. The numeric VS registration version contains only major and minor (currently 0.9), without the patch number. The Windows DLL product version is `0.9.0`, and its file version is `0.9.0.0`. KernelInfo retains its existing fields for computation diagnostics.

## Threads and execution parameters

Both filters finish each request in the host calling thread, creating no internal workers or thread pool. PocketFFT also uses one thread. SIMD and batched FFT remain enabled without requiring extra threads. The host controls concurrent requests for different output frames.

| Interface | Current behavior | Script migration |
|---|---|---|
| DFTTest threads, fft_threads | Optional integer compatibility arguments; never read, normalized or stored. | Existing arguments may remain; no workers are created. |
| FFT3D mt, ncpu, measure | Optional compatibility arguments: bool, int, bool; completely ignored. | Existing arguments may remain. |
| FFT3D / DFTTest fft_backend | Optional string compatibility argument; content ignored. PocketFFT remains fixed. | Existing arguments may remain; no backend name changes execution. |
| KernelInfo fft_backend output | Retained: actual FFT implementation and SIMD path. | Read-only diagnostic, not backend selection. |
| KernelInfo fft_threads output | Removed; FFT backend/profile/SIMD diagnostics remain. | Delete field reads; compatibility inputs add no output fields. |

Compatibility arguments exist only in host registration signatures and never enter DS2 parameter parsing or algorithm configuration. Hosts still enforce the registered types; the plugin does not validate their numeric ranges or string contents. VS compatibility integers are not restricted to int32. Newly restored arguments follow the current signature without moving existing parameters; check positional calls when migrating older signatures, or use keywords.

Each instance retains at most one idle workspace across planes, capped at 64 MiB. Active requests can acquire separate workspaces; this idle cap limits neither host concurrency nor total process memory. FFT3D raw spectra have separate cache_frames/cache_mb budgets. See [execution and precision](../../knowledge/en/shared/execution-precision.md).

## Load the plugin

```python
import vapoursynth as vs
core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
print(core.neo_fft.KernelInfo())
```

Replace the path: the Windows binary is neo-fft.dll, Linux neo-fft.so. Omit LoadPlugin in an autoload setup. Start the mathematical explanation with [FFT3D](../../knowledge/en/fft3d.md) or [DFTTest](../../knowledge/en/dfttest.md).

## AviSynth calls and builds

Requires an AviSynth+ host supporting interface version 11. Load the same neo-fft.dll with LoadPlugin, then call neo_fft_FFT3D, neo_fft_DFTTest or neo_fft_KernelInfo. Shared filter arguments retain VS semantics. AVS y/u/v/a plane modes sit between fft_backend and kalman_warmup in FFT3D, and after all shared arguments in DFTTest. Planar YUVA/RGBA is supported; Alpha is copied by default and can be selected explicitly. Arrays use native [0,1] syntax; a scalar is accepted as a one-element array. Both filters use `planes=[]` to select no planes and copy everything. Omitting the argument defaults to all existing non-Alpha planes; AVS can still select through y/u/v/a. To express omission explicitly, write `planes=None` on VS or `planes=Undefined()` on AVS (parentheses required). Omitted values do not override AVS y/u/v/a; `[]` does.

Both video filters forward input audio and parity; output properties still come from source frame n. KernelInfo returns a native array in the order documented on its API page. Filters remain internally single-threaded; Prefetch concurrency belongs to the host.

DFTTest sample locations and its four curve parameters also accept numeric strings separated by whitespace, commas or colons; see [DFTTest AviSynth calls](dfttest.md#avisynth).

NEO_FFT_BUILD_AVISYNTH and NEO_FFT_BUILD_VAPOURSYNTH default to ON. Disable either for a single-host build, or both for core only. Windows AVS C++ ABI requires MSVC or clang-cl; MinGW builds must disable AVS. Set NEO_FFT_AVS_SDK for an SDK override. Host tests additionally require NEO_FFT_TEST_AVISYNTH=ON and NEO_FFT_AVISYNTH_RUNTIME. Outputs are neo-fft.dll on Windows and neo-fft.so on Linux.
