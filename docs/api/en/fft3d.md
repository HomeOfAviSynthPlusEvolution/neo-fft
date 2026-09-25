# FFT3D

[Contents](README.md)

Window overlapping image blocks, transform them to frequency space, suppress noise and reconstruct the output. FFT3D also provides Kalman recurrence, sharpening and dehalo. The default uses three adjacent frames at the same spatial positions, without motion estimation or compensation.

## Calling convention

```text
core.neo_fft.FFT3D(
    clip [, sigma, beta, planes, bw, bh, bt, ow, oh,
    kratio, sharpen, scutoff, svr, smin, smax, interlaced,
    wintype, pframe, px, py, pshow, pcutoff, pfactor,
    sigma2, sigma3, sigma4, degrid, dehalo, hr, ht,
    l, t, r, b, opt, cache_frames, cache_mb, mt, ncpu, measure, fft_backend, kalman_warmup]
) -> VideoNode
```

This is registration order. Brackets denote optional arguments, not Python arrays. Only clip is required; keyword arguments are recommended. Input must have fixed dimensions and format, at least one frame, and planar GRAY/YUV/RGB with 1 or 3 planes (AviSynth also accepts four-plane YUVA/RGBA). Supported samples are 8/10/12/14/16-bit integers and float32. Every selected plane must pass the geometry checks below.

## Input, noise and grid compensation

| Argument | Type | Default | Domain and effect |
|---|---|---|---|
| `clip` | VideoNode | required | Input; dimensions, format, frame count and rate are preserved. |
| `planes` | integer array | all non-Alpha planes | Omitted or VS `None` selects all existing non-Alpha planes; explicit `[]` selects none and copies everything. Valid plane indices only; duplicates are processed once. Unselected planes are copied. |
| `sigma` | float | `2.0` | Nonnegative noise standard deviation in 8-bit units; multiplied by 256 for 16-bit, divided by 255 for float. |
| `beta` | float | `1.0` | At least 1. Wiener gain floor `(beta-1)/beta`; larger values retain more residual. Kalman does not use this gain but still validates the argument. |
| `sigma2`, `sigma3`, `sigma4` | float | each inherits `sigma` | Nonnegative frequency-profile controls together with sigma. Equal values select uniform noise; these are not per-temporal-frame noise estimates. |
| `degrid` | float | `1.0` | Nonnegative strength of removal/restoration of the window-shaped component estimated from block DC; 0 disables it. It does not merely remove DC. Kalman recurrence does not use degrid, but noise sampling and enabled enhancement can. |

Positive pfactor selects sampled noise instead of sigma1..4 for ordinary Wiener; sigma1 means sigma. Kalman initialization still uses sigma even with sampled noise.

## Blocks, time and processing region

| Argument | Type | Default | Domain and effect |
|---|---|---|---|
| `bw`, `bh` | integer | `32`, `32` | Each at least 2, in samples of each plane. Chroma uses the same block dimensions, not subsampling-scaled sizes. Width and height may differ. |
| `bt` | integer | `3` | Inclusive range -1..5; modes below. |
| `ow`, `oh` | integer | `-1`, `-1` | Any negative value derives floor(block length/3). Nonnegative overlap must be at most floor(block length/2). Default block 32 gives overlap 10, step 22. |
| `wintype` | integer | `0` | 0, 1 or 2: overlap shoulders use cosine/cosine, square-root cosine/its cube, or rectangular/cosine squared for analysis/synthesis. |
| `l`, `t`, `r`, `b` | integer | each `0` | Nonnegative full-resolution margins excluded from filtering. This defines an ROI, not output cropping. |
| `interlaced` | bool | `False` | Reorders ROI rows for field-style processing and restores their order. Each selected ROI height must be even. This is not deinterlacing and does not change frame rate. |

| bt | Mode | Ordinary temporal neighborhood for output n |
|---:|---|---|
| -1 | Frequency enhancement only | n |
| 0 | Kalman recurrence | Bounded warmup or nearby checkpoint continuation; frame 0 is copied |
| 1 | Spatial Wiener | n |
| 2 | Temporal Wiener | n−1, n |
| 3 | Temporal Wiener | n−1, n, n+1 |
| 4 | Temporal Wiener | n−2, n−1, n, n+1 |
| 5 | Temporal Wiener | n−2, n−1, n, n+1, n+2 |

For bt=2..5, an incomplete neighborhood makes the entire request use spatial T=1 processing. Endpoints are not replicated and no other shorter temporal length is substituted. Clips shorter than bt are legal. Sampled noise can additionally require pframe.

l+r must be less than frame width and t+b less than height. If YUV chroma is selected, margins must align to its subsampling on each axis; luma-only processing need not satisfy alignment for unselected chroma. Each ROI must contain a bw×bh block and satisfy the single-reflection padding limits. Merely fitting one block is insufficient: axis length 8, block 8 and overlap 0 fails padding. See [geometry](../../knowledge/en/shared/block-geometry.md).

## Noise sampling and preview

| Argument | Type | Default | Domain and effect |
|---|---|---|---|
| `pfactor` | float | `0.0` | Nonnegative; positive values select sampled noise and linearly scale sampled power for ordinary Wiener. For Kalman, positive values only select sampling and do not scale power. |
| `pframe` | integer | `0` | Noise sample frame, clamped to clip endpoints. Effective preview instead uses current output n. |
| `px`, `py` | integer | `0`, `0` | Nonnegative block-grid indices, not pixel coordinates. Both zero means automatic search; otherwise select a block manually. |
| `pcutoff` | float | `0.1` | Strictly positive high-pass scoring scale. |
| `pshow` | bool | `False` | Preview the selected noise block in the current frame when sampling or a nonuniform sigma profile is active. Ineffective for uniform noise. |

Priority is positive pfactor → unequal effective sigma values → uniform sigma. Automatic selection excludes two block layers on every edge: bx=2..Nx−3, by=2..Ny−3. Every relevant grid must therefore be at least 5×5. Manual indices must fit each relevant plane; chroma grids are usually smaller.

Effective preview overrides bt, denoising, Kalman and enhancement. It reads neither temporal neighbors nor pframe and does not advance Kalman state. It shows the selected block's intersection with the visible ROI; the rest of the ROI has zero reconstruction, with the integer chroma midpoint restored. Outside selected ROIs and on unselected planes input is copied, except explicit AVS mode 1. No text overlay is drawn.

## Kalman and enhancement

| Argument | Type | Default | Domain and effect |
|---|---|---|---|
| `kratio` | float | `2.0` | Nonnegative spectral-change threshold factor for Kalman resets. Validated even in other modes. |
| `kalman_warmup` | integer | `8` | 0..2147483647; maximum historical recurrence frames excluding target n. Only effective bt=0 uses it; all modes validate it. 0 disables cold pre-roll, while sequential requests can still reuse state. No auto or unlimited sentinel. |
| `sharpen` | float | `0.0` | Nonnegative sharpening strength; 0 disables. |
| `scutoff` | float | `0.3` | Positive transition scale for sharpening frequency weights. |
| `svr` | float | `1.0` | Nonnegative vertical-frequency scale for enhancement; 0 is valid. |
| `smin`, `smax` | float | `4.0`, `20.0` | 0≤smin≤smax; scaled like sigma before computing sharpening gain. smax=0 gives no sharpening gain. |
| `dehalo` | float | `0.0` | Nonnegative dehalo strength; 0 disables. |
| `hr` | float | `2.0` | Positive dehalo frequency scale. With dehalo active, sampled weights must also have a finite positive normalization maximum. |
| `ht` | float | `50.0` | Nonnegative dehalo power scale. Unlike sigma and smin/smax, it is not bit-depth-scaled. |

Enhancement generally operates on recovered current-frame spatial spectra. Uniform single-frame Wiener fuses its gain with enhancement using a specific shared power calculation; stages cannot be arbitrarily reordered. bt=-1 with both enhancement strengths zero still performs windowing, FFT and reconstruction, rather than a byte copy.

### Bounded Kalman warmup

AVS and VS both accept `kalman_warmup=8`; omission, VS `None` and AVS `Undefined()` use the default. The argument is appended after the old host slots (after y/u/v/a on AVS).

For n>0, reuse the nearest completed checkpoint j≤n only if n-j≤W+1. Otherwise initialize L=0,C=Q=R0 and consume max(1,n-W)..n. The default cold request for frame 10000 processes 9992..10000 (9 frames); checkpoint 9996 reduces this to 9997..10000 (4); checkpoint 9999 reduces it to just 10000. Sequential requests normally reuse long history instead of resetting every W frames. Frame 0 remains a copy except explicit AVS mode 1.

The limit covers recurrence work. A cold sampled-noise model may additionally fetch pframe outside the window, so W+1 is not a total source-callback or wall-clock guarantee. Spectral motion resets are not whole-frame scene detection.

This intentionally trades request-order independence for bounded cold work. Seeks, eviction and concurrent scheduling can choose different histories and produce different results; even repeated frames can differ after eviction. Private state prevents concurrent mutation, not numerical history differences. Short warmup may temporarily reduce denoising; increasing it costs more upstream work and does not guarantee convergence. See [Kalman state and checkpoints](../../knowledge/en/fft3d/kalman.md).

```python
out = core.neo_fft.FFT3D(clip, bt=0, kalman_warmup=8)
```

```avs
out = neo_fft_FFT3D(clip, bt=0, kalman_warmup=8)
```

## Execution and spectrum cache

| Argument | Type | Default | Domain and effect |
|---|---|---|---|
| `opt` | integer | `0` | 1 selects this project's scalar kernels; other int32 values select automatic SIMD, or scalar if SIMD is not built. Values do not force a specific ISA. opt=1 does not force the FFT library to scalar. |
| `cache_frames` | integer | `-1` | -1 derives bt + host concurrency − 1: VS reads core threads at creation; AVS uses thread cache hints, assuming 1 until notified. 0 disables; positive values limit distinct source frame IDs retained, not preallocated full frames. |
| `cache_mb` | integer | `128` | Positive MiB budget, 1 MiB=1048576 bytes. -1 also uses 128; 0 disables. Limits charged cache entries, not process memory. |

FFT3D creates no internal workers: selected planes run sequentially in the host calling thread. VS may concurrently request output frames. FFT is fixed to PocketFFT with one internal thread. Optional compatibility arguments mt (bool), ncpu (int), measure (bool) and fft_backend (string) are accepted by the host but never read, value-validated or stored by the plugin. They cannot change threading, planning or the backend; see [compatibility arguments](README.md#threads-and-execution-parameters).

Both cache limits apply; either 0 disables caching, values below -1 fail. There is no unlimited sentinel. Larger positive values remain finite limits. Only ordinary processing configured with bt=2..5 uses this cache, including its spatial endpoint fallback. Effective preview and bt=-1/0/1 do not use it.

Insufficient budget falls back to private block computation without waiting for capacity or changing the filter mode. Concurrent requests can share an admitted row. These arguments do not control Kalman checkpoints, VS frame caching or active workspaces. See [spectrum caching](../../knowledge/en/fft3d/spectra-cache.md).

## Return value and frame properties

Returns a VideoNode preserving dimensions, format, length and rate; output properties come from current input frame n. Unselected planes (except explicit AVS mode 1) and areas outside selected ROIs are copied byte for byte. No cache-statistics properties are added.

Integer conversion restores chroma midpoint, rounds and clips to the bit-depth range. **Processed float ROI output is clipped to [0,1], including float YUV chroma**, so negative chroma can become zero. Copied regions preserve their values. Kalman output frame 0 copies every plane except explicit AVS mode 1 and skips these conversions.

## Minimal examples

Replace the plugin path; BlankClip needs no source plugin.

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
clip = core.std.BlankClip(width=128, height=96, length=12,
                          format=vs.YUV420P8, color=[96, 128, 128])
output = core.neo_fft.FFT3D(clip, sigma=2.0, bt=3)
output.set_output()
```

The alternatives below reuse core and clip from above:

```python
# Luma only, five-frame Wiener, at most 256 MiB of row spectra.
output = core.neo_fft.FFT3D(clip, planes=[0], bt=5, cache_mb=256)
output.set_output()
```

```python
# Kalman recurrence; cache_mb does not control its checkpoints.
output = core.neo_fft.FFT3D(clip, bt=0, sigma=2.0, kratio=2.0)
output.set_output()
```

## Limits and common errors

| Situation | Behavior or remedy |
|---|---|
| Explicitly request default plane selection | Use `planes=None` on VS or `planes=Undefined()` on AVS; `planes=[]` copies everything. |
| Large blocks on small chroma/ROI | Block or reflection geometry may fail; reduce blocks or select suitable planes. |
| Pass pixel coordinates as px/py | Use block-grid indices derived separately for each selected plane. |
| Missing temporal neighbors / clip shorter than bt | Valid bt=2..5 falls back to spatial processing; invalid bt still fails. |
| Invalid bt, wintype, beta, planes; negative sigma/margins | Creation fails even if the current mode does not consume the argument. |
| NaN, Inf or numeric parameters outside float32 range | Creation fails. Integer arguments and array elements must fit int32, except completely ignored execution compatibility arguments. |
| Finite but extreme parameters or pixels | Active arithmetic can overflow and fail creation/frame processing; partial output is not published. |
| Consumed float ROI contains NaN/Inf | Frame processing fails. Pure copies are not scanned, including Kalman frame 0. |
| Treat interlaced as deinterlacing or field-order detection | It only reorders ROI rows; it neither infers TFF/BFF nor rewrites _FieldBased. |
| Treat cache budget as peak process memory | Active workspaces, source frames, host caches and other state require additional memory. |

See [the processing pipeline](../../knowledge/en/fft3d.md), [noise and enhancement](../../knowledge/en/fft3d/noise-enhancement.md), [Kalman](../../knowledge/en/fft3d/kalman.md) and [execution policy](README.md#threads-and-execution-parameters).

## Omitted arguments and empty lists

On VS, explicit `planes=None` is equivalent to omitting the argument. On AVS, use `planes=Undefined()` with parentheses. Neither is an empty list: in both filters, `planes=[]` selects no planes and copies the current frame bit for bit after configuration validation. The default processes all existing non-Alpha planes. On AVS, omission or `Undefined()` still uses `y/u/v/a` and their defaults. An explicit empty list overrides `y/u/v/a`, including mode 1 or invalid mode values.

```python
core.neo_fft.FFT3D(clip, planes=None)  # Default selection
core.neo_fft.FFT3D(clip, planes=[])    # Copy everything
```

```avs
neo_fft_FFT3D(c, planes=Undefined(), y=3, u=2, v=2) # Use y/u/v/a
neo_fft_FFT3D(c, planes=[])                         # Copy everything
```

## AviSynth

Both filters accept planar YUVA/RGBA, with Alpha at index 3. Alpha is full resolution and uses luma/RGB sample scaling, not chroma centering or subsampling. Omitted selection processes existing non-Alpha planes and copies Alpha. Explicit `planes=[3]` processes Alpha alone.

The AVS-only integer arguments `y`, `u`, `v`, `a` sit in that order between `fft_backend` and `kalman_warmup`. If `planes` is supplied, it takes precedence and these modes are ignored, including when `planes=[]`. Otherwise omitted modes default to `y=u=v=3, a=2`; mode 1 skips writes, mode 2 copies, and mode 3 processes. Mode 1 leaves the output plane uninitialized with unspecified pixel contents, for workflows that later extract and combine only useful planes. It avoids copying that plane but does not reduce the output format or allocation. Valid modes for missing planes have no effect; other mode values fail unless overridden by `planes`. On RGB these names select R/G/B/A respectively. `y=2,u=2,v=2,a=2` copies every plane, including for FFT3D. VS retains its `planes` interface; Alpha held as a separate gray clip can be filtered separately.


AviSynth keeps the shared parameter names, inserts y/u/v/a at the position described above, and forwards input audio and parity. See [host requirements](README.md#avisynth-calls-and-builds).

```avs
LoadPlugin("/path/to/neo-fft.dll")
c=BlankClip(width=128,height=96,length=12,pixel_type="YV12",color_yuv=$608080)
return neo_fft_FFT3D(c,sigma=2.0, bt=3)
```
