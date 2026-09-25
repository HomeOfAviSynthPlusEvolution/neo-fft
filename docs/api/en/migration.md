# Migrating from older Neo filters

[Contents](README.md)

This guide covers AviSynth / VapourSynth scripts using older **Neo FFT3D** (neo-fft3d) and **Neo DFTTest** (neo-dfttest). neo-fft combines them into one plugin with separate algorithm entry points. It does not promise compatibility with every historical interface of traditional FFT3DFilter or tritical/pinterf DFTTest; for example, nfile, nstring and sstring are outside this same-name parameter migration.

Common CPU functionality can be migrated, but accepting a call and reproducing the old output pixel for pixel are separate conditions. Required edits come first, followed by behaviors that can change results. Full parameter references: [FFT3D](fft3d.md), [DFTTest](dfttest.md) and [KernelInfo](kernel-info.md).

## Required edits

Apply only the rows relevant to your script. Explicit loading paths, function names or namespaces must use the new entry points.

| Situation | Old call | New call or required action |
|---|---|---|
| Explicit plugin loading | neo-fft3d.dll, neo-dfttest.dll | Load neo-fft.dll; the Linux binary is neo-fft.so. Install the new plugin when using autoload. |
| AVS FFT3D | `neo_fft3d(c, ...)` | `neo_fft_FFT3D(c, ...)` |
| AVS DFTTest | `neo_dfttest(c, ...)` | `neo_fft_DFTTest(c, ...)` |
| VS FFT3D | `core.neo_fft3d.FFT3D(c, ...)` | `core.neo_fft.FFT3D(c, ...)` |
| VS DFTTest | `core.neo_dfttest.DFTTest(c, ...)` | `core.neo_fft.DFTTest(c, ...)` |
| Tools looking up plugins by identifier | Old plugin identifier | Use `org.neofilters.neo_fft`. |
| Positional arguments after clip | Depend on the old signature order | Rewrite according to the new API, preferably using named arguments. Slots are not guaranteed to match the old signature. |
| AVS DFTTest on YUVA/RGBA, preserving the old default Alpha processing | Omitted `a`, formerly defaulting to 3 | Add `a=3` when `planes` is absent; otherwise include index 3 in the explicit `planes` selection. |
| Old diagnostic fields or return structures | Reading thread counts or other old fields | Adapt to KernelInfo: VS returns a dictionary; AVS returns a four-element array. Thread count fields are no longer returned. |

AVS requires an AviSynth+ host supporting interface version 11; older hosts need upgrading. Windows AVS builds use MSVC or clang-cl. See the [API index](README.md) for host build and loading requirements.

These examples assume an existing input `c`. They replace the old entry point and use named arguments:

```avs
LoadPlugin("/path/to/neo-fft.dll")
# Old call: neo_fft3d(c, sigma=2.0, bt=3, y=3, u=2, v=2)
return neo_fft_FFT3D(c, sigma=2.0, bt=3, y=3, u=2, v=2)
```

```python
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
# Old call: core.neo_dfttest.DFTTest(c, sigma=8.0, tbsize=3, planes=[0])
output = core.neo_fft.DFTTest(c, sigma=8.0, tbsize=3, planes=[0])
```

### Arguments you can keep, and controls that no longer take effect

AVS `y/u/v/a` and VS `planes` remain available; AVS also accepts `planes`. Explicit `planes` overrides every legacy plane mode, even when it is an empty array; both filters copy everything for `planes=[]`. VS `planes=None` and AVS `planes=Undefined()` (parentheses required) mean omission and do not override AVS y/u/v/a. With no selection arguments, both filters process existing non-Alpha planes and copy Alpha. Without `planes`, omitted modes default individually to `y=u=v=3, a=2`; 1 skips writes, 2 copies, and 3 processes. Mode 1 matches the old interface: output pixels in that plane are unspecified, allowing scripts to continue extracting and combining only useful planes. On RGB these names map to R/G/B/A.

The verified old Neo FFT3D r14 and Neo DFTTest r12 VS interfaces both reject explicit `planes=[]`. Supporting an empty selection is an extension, not a required edit for those old scripts. Earlier neo-fft development builds treated FFT3D `[]` as the default selection; scripts relying on that behavior should omit the argument or use VS `None` / AVS `Undefined()`.

AVS DFTTest accepts strings separated by whitespace, commas, colons or mixtures of them for `nlocation`, `slocation`, `ssx`, `ssy` and `sst`, as well as native numeric arrays. Existing strings do not need delimiter changes. VS continues to take numeric arrays. The decimal point is always `.`, and commas always separate values. For example, `slocation="0:4, 1:12"` is equivalent to `slocation=[0,4,1,12]`.

FFT3D `mt/ncpu/measure`, DFTTest `threads/fft_threads`, and both filters' `fft_backend` may remain. The host checks their registered types, but the plugin never reads their contents. They create no internal workers and do not switch the fixed PocketFFT backend. Scripts relying on these arguments to allocate parallel work must instead configure host frame concurrency.

`opt=1` selects this project's scalar kernels; other admitted values use automatic SIMD. These values no longer select or restrict a particular ISA, and `opt=1` does not force the FFT implementation to scalar. FFT3D accepts any int32; DFTTest accepts 0/1/2/3/8. Scripts relying on old opt values as ISA restrictions cannot treat them as equivalent controls; see each filter's API.

## Differences in results

These include intentional behavior definitions and numerical exceptions found in reference comparisons. They cannot all be classified as ordinary FFT rounding. Reference differences concern the tested old builds, not necessarily every historical release.

| Situation | Current behavior and migration effect |
|---|---|
| Default AVS Alpha selection | Alpha is copied; older Neo DFTTest processed it by default. Select it explicitly to preserve processing, as described above. Alpha uses full resolution and the luma/RGB sample domain. |
| Positive FFT3D `pfactor` in ordinary Wiener filtering | Linearly scales sampled noise power; the fixed old reference only distinguished zero from nonzero. Positive values other than 1 may change denoising strength. Use `pfactor=1` on both sides for matched comparisons; other strengths remain useful. Kalman is an exception: positive values select sampling without scaling its power. |
| FFT3D `pshow` | Uses separate preview windows and the actual plane's chroma flag, avoiding the old reference's window-state and luma-offset issues. Repeated or reordered previews and some plane outputs may differ; those old behaviors are not reproduced. |
| FFT3D Kalman | Uses bounded warmup (kalman_warmup=8) and nearby completed checkpoints. Cold seeks do not replay the full prefix, and results can depend on cache/request history or concurrency; this differs from the old implementation's call-history policy and the former neo-fft canonical replay. Even sequential requests can differ due to rounding near reset thresholds; see below. |
| DFTTest `dither>=2` | Noise depends on seed, frame, plane and pixel coordinates rather than following the old implementation's evolving worker-slot or request-order sequence. The same seed does not guarantee the old noise pattern. |
| Floating-point arithmetic and quantization | FFT implementation, SIMD and operation grouping can change rounding. Integer output can also differ at thresholds or final quantization. Historical test tolerances cover their specific inputs and configurations, not every possible input. |
| Invalid or degenerate configurations | Non-finite values, duplicate curve knots and out-of-bounds samples are explicitly rejected. Algorithm parameters may still require raw validation when inactive in the current mode. Invalid indexing and division-by-zero behavior are not reproduced; explicit plane mode 1 retains its skip-write semantics. Completely ignored execution compatibility arguments are exempt from these algorithm-value checks. |

### Known Kalman reference exception

In the existing AVS comparison against fixed reference `bc9dd394`, six configurations (B32, sigma=2, kratio=2, 128×96 RGB16/float, covering three dispatch modes) exceeded the frozen pixel tolerances: RGB16 maximum difference 33 against a budget of 4; float maximum about 7.9751e-4 against 2e-5. They remain reference mismatches. The tolerances were not widened, and these cases must not be counted as passes.

The divergence was traced to windowing/FFT rounding crossing a strict motion threshold, followed by amplification through recurrence and reconstruction. No violation of the defined operators was found in the current implementation. This does not mean every Kalman input has that error, or that every old Kalman implementation is mathematically wrong. Users who need to preserve existing pixel-exact golden outputs should check their Kalman configurations separately.

Background in the same-language knowledge base: [noise and enhancement](../../knowledge/en/fft3d/noise-enhancement.md), [Kalman state and replay](../../knowledge/en/fft3d/kalman.md), [dither](../../knowledge/en/dfttest/dither.md) and [execution and precision](../../knowledge/en/shared/execution-precision.md).
