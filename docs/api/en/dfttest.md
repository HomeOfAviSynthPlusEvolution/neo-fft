# DFTTest

[Contents](README.md)

Filters colocated space/time blocks in the frequency domain, with five spectral gains, curves, sampled noise models and temporal overlap-add. It does not estimate motion vectors.

## Invocation

```text
core.neo_fft.DFTTest(
    clip [, ftype, sigma, sigma2, pmin, pmax, sbsize, smode, sosize,
    tbsize, tmode, tosize, swin, twin, sbeta, tbeta, zmean, f0beta,
    nlocation, alpha, slocation, ssx, ssy, sst, ssystem, dither,
    dither_seed, planes, opt, threads, fft_threads, fft_backend]
) -> VideoNode
```

Brackets mark optional arguments; only clip is required. This is registration order; prefer keyword arguments. Input must have fixed dimensions/format and at least one frame: planar GRAY/YUV/RGB, one or three planes (four-plane YUVA/RGBA on AviSynth), integer 8/10/12/14/16-bit or float32. The AviSynth function is neo_fft_DFTTest, with the same shared argument order, followed by AVS-only y/u/v/a.

## Filtering controls

| Parameter | Type | Default | Domain and effect |
|---|---|---|---|
| `clip` | VideoNode | required | Input clip |
| `ftype` | int | 0 | 0 Wiener; 1 hard threshold; 2 direct gain; 3 separate gains inside/outside a power interval; 4 smooth power-dependent gain |
| `sigma` | float | 8 | Nonnegative; a power parameter for 0/1, a gain for 2/3/4. **Not squared**; FFT3D's standard-deviation interpretation does not apply |
| `sigma2` | float | 8 | Nonnegative; outside-interval gain for type 3. Independently defaults to 8, not sigma |
| `pmin`, `pmax` | float | 0, 500 | `0≤pmin≤pmax`; power parameters for types 3/4, calibrated by window energy |
| `f0beta` | float | 1 | Positive exponent for the type-0 gain; unused by other types |
| `zmean` | bool | True | Remove and restore a DC-derived window-shaped spectrum over the whole space/time volume, not per frame |

## Blocks and windows

| Parameter | Type | Default | Domain and effect |
|---|---|---|---|
| `sbsize` | int | 16 | Positive, in samples of each plane; must be odd for smode=0 |
| `smode` | int | 1 | 0 outputs only the center sample of each pixel-centered block; 1 spatial overlap-add |
| `sosize` | int | 12 | smode=0 normalizes any int32 to 0. Otherwise `0≤O<S`; if `O>floor(S/2)`, require `S%(S-O)=0`, S=sbsize |
| `tbsize` | int | 3 | `1..15`, no greater than clip length; odd for tmode=0, either parity for tmode=1 |
| `tmode` | int | 0 | 0 centered window and center-only temporal output; 1 fixed temporal lattice, combining every block covering the output |
| `tosize` | int | 0 | tmode=0 normalizes any int32 to 0. Otherwise `0≤O<T`; if `O>floor(T/2)`, require `T%(T-O)=0` |
| `swin`, `twin` | int | 0, 7 | Spatial/temporal raw windows, 0..11; evaluated even at length 1 |
| `sbeta`, `tbeta` | float | 2.5, 2.5 | Nonnegative Kaiser parameters; unused by other corresponding windows |

Window IDs: 0 Hann, 1 Hamming, 2 Blackman, 3 four-term Blackman–Harris, 4 Kaiser, 5 seven-term Blackman–Harris, 6 flat top, 7 rectangular, 8 Bartlett, 9 asymmetric Bartlett–Hann variant, 10 Nuttall, 11 Blackman–Nuttall. Names do not replace the exact [window equations](../../knowledge/en/shared/windows-reconstruction.md). Active windows need finite nonzero energy; zmean additionally needs nonzero template DC.

Spatial padding uses one endpoint-excluding reflection, checked separately for each selected plane, including small chroma planes. Temporal endpoints repeat the nearest frame while retaining T logical slots; T must still fit clip length. `tmode=1,T=4,O=2` is valid; `T=5,O=3` is not. Omitted tosize always remains 0.

## Curves and sampling

| Parameter | Type | Default | Domain and effect |
|---|---|---|---|
| `slocation` | float[] | empty | Shared frequency curve; overrides axis curves when nonempty |
| `ssx`, `ssy`, `sst` | float[] | empty each | X, Y and temporal curves; empty axes fall back to sigma |
| `ssystem` | int | 0 | 0 separable axis product; 1 radial. Prefer slocation for radial profiles; with axis-only curves the temporal table supplies radial ordinates |
| `nlocation` | int[] | empty | Up to 500 flattened `(start_frame, plane, y, x)` tuples, each selecting an S×S raw patch over T consecutive frames |
| `alpha` | float | 5 for ftype=0, otherwise 7 | Positive sampled-power multiplier |

Curves are `[frequency,value,...]`. Each nonempty array needs endpoints 0 and 1, positions in `[0,1]`, nonnegative values, and positions unique after float32 conversion. Input may be unsorted. Overridden curves still undergo structural validation. Values are transformed before interpolation according to the selected system; see [spectral models](../../knowledge/en/dfttest/spectral-models.md).

Sampling requires `0≤fn≤N-T`, an existing plane, `0≤x≤plane_width-S` and `0≤y≤plane_height-S`. Coordinates refer to actual plane pixels, not luma coordinates or block indices; sampling does not reflect or clamp. Duplicate tuples retain their averaging weight; a tuple may select an output-unselected plane. For types 0/1, nonempty nlocation **replaces** the primary sigma/curve table. Types 2/3/4 still validate sampling parameters but request no sampling frames. All selected output planes share one sampled table.

## Output and execution

| Parameter | Type | Default | Domain and effect |
|---|---|---|---|
| `dither` | int | 0 | Nonnegative; selected UInt8 planes only. 0 ordinary conversion; 1 error diffusion; ≥2 diffusion plus deterministic noise |
| `dither_seed` | int | 0 | Nonnegative; affects noise only for UInt8 with dither≥2; unused by 0/1 |
| `planes` | int[] | all non-Alpha when omitted | Explicit `[]` selects none; valid actual indices, duplicates selected once |
| `opt` | int | 0 | Accepts 0/1/2/3/8; 1 selects own scalar kernels, others automatic Highway. Does not force an ISA or scalar FFT |
| `threads`, `fft_threads` | int | omitted | Call compatibility only; never read, normalized or stored, with no plugin int32 restriction or execution effect |
| `fft_backend` | string | omitted | Call compatibility only; content ignored, no backend selection |

Except for the ignored compatibility arguments above, all integers/array integers must fit int32, and floating parameters must be finite and representable as float32. Inactive algorithm parameters are still validated. Compatibility arguments are subject only to host signature types and never enter DS2 parsing or algorithm configuration. PocketFFT is fixed internally; neither filter nor FFT creates workers. The host can process concurrent frame requests. There is no DFTTest cross-frame spectrum-cache parameter.

The returned clip preserves dimensions, format, frame count and rate; properties come from source frame n. Unselected planes are bitwise copies, except explicit AVS mode 1. After configuration validation, `planes=[]` copies n without constructing models or FFT windows. Integer output is rounded and clipped once. **Float output is divided by 255 without clipping**, retaining negative chroma and out-of-nominal values. Dither runs once after all space/time contributions; other bit depths ignore valid dither settings.

## Minimal examples

```python
import vapoursynth as vs
core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
clip = core.std.BlankClip(width=128, height=96, length=12,
                          format=vs.YUV420P8, color=[96, 128, 128])
output = core.neo_fft.DFTTest(clip, sigma=8.0, tbsize=3, planes=[0])
output.set_output()
```

Each following call can replace output using that source:

```python
# Four-frame volumes with two-frame temporal overlap.
output = core.neo_fft.DFTTest(clip, tmode=1, tbsize=4, tosize=2)
# Shared radial frequency curve, default spatial block/overlap.
output = core.neo_fft.DFTTest(clip, slocation=[0, 4, 1, 12], ssystem=1)
# Upper-left luma 16x16 patch over three frames starting at frame 0.
output = core.neo_fft.DFTTest(clip, nlocation=[0, 0, 0, 0], alpha=5)
```

## Common errors

Reduce sosize when reducing sbsize: smode=1,sbsize=8 cannot retain the default overlap 12. Center mode may use sbsize=5 without changing sosize. T exceeding clip length, even T in tmode=0, invalid overlap, zero-energy windows or zero-DC mean templates fail creation. Consumed float samples or intermediate values containing NaN/Inf fail the frame; pure copies of unselected regions are not scanned for finiteness.

See [DFTTest computation](../../knowledge/en/dfttest.md), [temporal overlap](../../knowledge/en/dfttest/temporal-ola.md) and [dither](../../knowledge/en/dfttest/dither.md).

## Omitted arguments and empty lists

On VS, explicit `planes=None` is equivalent to omitting the argument. On AVS, use `planes=Undefined()` with parentheses. Neither is an empty list: in both filters, `planes=[]` selects no planes and copies the current frame bit for bit after configuration validation. The default processes all existing non-Alpha planes. On AVS, omission or `Undefined()` still uses `y/u/v/a` and their defaults. An explicit empty list overrides `y/u/v/a`, including mode 1 or invalid mode values.

```python
core.neo_fft.DFTTest(clip, planes=None)  # Default selection
core.neo_fft.DFTTest(clip, planes=[])    # Copy everything
```

```avs
neo_fft_DFTTest(c, planes=Undefined(), y=3, u=2, v=2) # Use y/u/v/a
neo_fft_DFTTest(c, planes=[])                         # Copy everything
```

## AviSynth

Both filters accept planar YUVA/RGBA, with Alpha at index 3. Alpha is full resolution and uses luma/RGB sample scaling, not chroma centering or subsampling. Omitted selection processes existing non-Alpha planes and copies Alpha. Explicit `planes=[3]` processes Alpha alone.

The AVS-only integer arguments `y`, `u`, `v`, `a` are appended in that order after all shared arguments. If `planes` is supplied, it takes precedence and these modes are ignored, including when `planes=[]`. Otherwise omitted modes default to `y=u=v=3, a=2`; mode 1 skips writes, mode 2 copies, and mode 3 processes. Mode 1 leaves the output plane uninitialized with unspecified pixel contents, for workflows that later extract and combine only useful planes. It avoids copying that plane but does not reduce the output format or allocation. Valid modes for missing planes have no effect; other mode values fail unless overridden by `planes`. On RGB these names select R/G/B/A respectively. `y=2,u=2,v=2,a=2` copies every plane, including for FFT3D. VS retains its `planes` interface; Alpha held as a separate gray clip can be filtered separately.


AviSynth keeps the shared parameter names and order, appending y/u/v/a after the execution compatibility arguments, and forwards input audio and parity. See [host requirements](README.md#avisynth-calls-and-builds).

`nlocation`, `slocation`, `ssx`, `ssy` and `sst` also accept numeric strings separated by whitespace, commas or colons and can be mixed with native arrays. `nlocation` takes decimal int32 integers; curves take decimal floating-point numbers, including signs and scientific notation. The decimal point is always `.`, independent of the system locale. Separators may be mixed, repeated, leading or trailing; empty strings or strings containing only separators mean empty arrays. A comma separates numbers, never a decimal fraction: "1,5" means two values, 1 and 5. Brackets, non-finite values and trailing junk are rejected. The same tuple, curve endpoint, coordinate and value checks described above apply after parsing. `planes` still takes integers or arrays; VS continues to use numeric arrays for these parameters.

```avs
LoadPlugin("/path/to/neo-fft.dll")
c=BlankClip(width=128,height=96,length=12,pixel_type="YV12",color_yuv=$608080)
return neo_fft_DFTTest(c,sigma=8.0, tbsize=3, planes=[0])
```

```avs
# Reuse c above; equivalent to slocation=[0,4,1,12].
return neo_fft_DFTTest(c,slocation="0:4, 1:1.2e1",ssystem=1)
# Sampling can likewise use nlocation="0,0:0,0".
```
