# DFTTest: phase-1 VS interface

Historical phase-1 contract. For current thread controls and removed parameters (including the fft_backend input), apply the later [execution-policy revision](../../phase-4/execution.md).

Specification: DFT-VS-001. This document defines the phase-1 registration, parameter defaults, supported subset, validation and output behavior.

## Interface

`core.neo_fft.DFTTest(clip, ftype, sigma, sigma2, pmin, pmax, sbsize, smode, sosize, tbsize, tmode, tosize, swin, twin, sbeta, tbeta, zmean, f0beta, nlocation, alpha, slocation, ssx, ssy, sst, ssystem, dither, dither_seed, planes, opt, threads, fft_threads, fft_backend) -> clip`

Arguments are in VS registration order; AVS-only y/u/v/a are absent. Only clip is signature-required; other arguments use the defaults below. It shares the `in.7086.neo_fft` identity and `neo_fft` namespace with FFT3D. Integer/array values must fit int32; floating values must remain finite through double-to-float conversion. Validate types before use, and validate finite values even on inactive parameters.

| Parameter | Type/default | Phase-1 domain and behavior |
| --- | --- | --- |
| clip | Required node | Fixed-size/format GRAY/YUV/RGB, 1 or 3 planes, UInt8/10/12/14/16 or Float32; positive dimensions and frame count |
| ftype | int 0 | 0..4; all implemented |
| sigma | float 8 | >=0; type-dependent threshold/gain, not standard deviation squared |
| sigma2 | float 8 | >=0; independently defaults to 8, **does not inherit sigma** |
| pmin, pmax | float 0,500 | 0<=pmin<=pmax; type 3/4 power parameters |
| sbsize | int 16 | >=1; odd required for smode=0 |
| smode | int 1 | 0=center sample; 1=overlap-add |
| sosize | int 12 | smode=0 forces effective 0; otherwise 0..sbsize-1 plus heavy-overlap divisibility rule |
| tbsize | int 3 | Only 1 implemented; omission resolves to 3 and fails explicitly |
| tmode | int 0 | Must be 0; reference also rejects 1 |
| tosize | int 0 | Effective 0 for tmode=0, regardless of supplied int32 value, following reference normalization |
| swin | int 0 | 0..11, exact equations in window spec |
| twin | int 7 | 0..11; evaluated at length 1, not silently ignored |
| sbeta, tbeta | float 2.5,2.5 | >=0; Kaiser parameters; inactive for other corresponding windows |
| zmean | bool true | Window-shaped mean removal/restoration |
| f0beta | float 1 | >0; exponent for ftype=0, inactive otherwise |
| nlocation | int array, empty | Must be empty; noise estimation deferred |
| alpha | float 5 if ftype=0, else 7 | Require that default while noise estimation is deferred |
| slocation, ssx, ssy, sst | float arrays, empty | Each must be empty; noise curves deferred |
| ssystem | int 0 | Must be 0 in this batch |
| dither | int 0 | Must be 0 |
| dither_seed | int, omitted | Omitted or nonnegative; accepted but inactive when dither=0, never initializes a hidden RNG |
| planes | int array, omitted | Omitted selects all actual planes; **explicit empty selects none**; indices 0..plane_count-1; duplicates select once |
| opt | int 0 | 0 automatic Highway; 1 own scalar kernels; other values deferred |
| threads | int 0 | <=0 automatic ->1; 1 explicit; >1 deferred |
| fft_threads | int 0 | <=0 automatic ->1; 1 explicit; >1 deferred |
| fft_backend | string `pocketfft` | RUN-001 backend policy; historical default was fftw |

`smode=0,sbsize=5` is valid without changing the default sosize=12 because normalization to zero precedes overlap validation. For smode=1,sbsize=8 the default sosize=12 is invalid and must be explicitly reduced. An unsupported feature request still fails when `planes=[]`; selecting no planes does not make an invalid configuration valid.

## Frame computation

1. At creation parse/normalize parameters, select planes and validate [geometry](kernel-geometry.md). Build [window/calibration](kernel-window-reconstruction.md) tables and [filter](kernel-filter.md) constants for selected planes.
2. Request only clip[n]. Preserve output metadata and copy unselected planes.
3. Reflect each selected plane, gather blocks in Y-then-X order and convert to windowed 8-bit amplitude units.
4. Forward FFT, optionally remove mean, run ftype, restore mean, then normalized inverse FFT.
5. Apply the explicit N compensation and synthesis; overlap-add or emit centers, crop and convert output without dither.

No temporal or noise-sample dependency is requested. Empty plane selection makes a bitwise visible-plane copy with source properties after configuration validation. Output size/format/count/rate do not change. Concurrency and failures follow [RUN-001](../execution.md).

## Examples and errors

`core.neo_fft.DFTTest(src, tbsize=1, fft_backend="pocketfft", opt=1)` retains ftype=0,sigma=8,sbsize=16,sosize=12,zmean=true. For GRAY8 128x96, P=152,Q=120,offset=(12,12),step=4, giving 35x27 block origins.

`core.neo_fft.DFTTest(src, tbsize=1, smode=0, sbsize=5, ftype=2, sigma=1, zmean=False, swin=7)` is a useful center-sample identity case, subject to admitted geometry and float rounding. `planes=[]` copies all planes, unlike FFT3D's empty-list behavior.

Reject even sbsize with smode=0, invalid overlaps/windows, unsafe padding, invalid planes, non-finite values and overflow at creation. Reject tbsize!=1, nonempty noise arrays, dither!=0 or unsupported execution controls as explicit unsupported requests. Non-finite selected frame content fails evaluation. Match errors by stage/category, not historical message spelling or old undefined behavior.

Acceptance is [P1-ACCEPT-001](../acceptance.md); all five ftypes, both modes and all spatial windows belong to this batch, not optional follow-up work.
