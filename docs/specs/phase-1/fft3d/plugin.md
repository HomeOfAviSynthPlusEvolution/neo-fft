# FFT3D: phase-1 VS interface

Historical phase-1 contract. For current thread controls and removed parameters (including the fft_backend input), apply the later [execution-policy revision](../../phase-4/execution.md).

Specification: F3D-VS-001. This document defines the phase-1 registration, parameter defaults, supported subset, validation and output behavior.

## Interface

`core.neo_fft.FFT3D(clip, sigma, beta, fft_backend, planes, bw, bh, bt, ow, oh, kratio, sharpen, scutoff, svr, smin, smax, measure, interlaced, wintype, pframe, px, py, pshow, pcutoff, pfactor, sigma2, sigma3, sigma4, degrid, dehalo, hr, ht, l, t, r, b, opt, ncpu, mt) -> clip`

Arguments above are in registration order. Only clip is required by the signature. All other arguments are optional with defaults below. The plugin identity is `in.7086.neo_fft`, namespace `neo_fft`; it also registers DFTTest. No AVS entry is required in this batch. Use named arguments in differential scripts.

Integers must fit int32 without narrowing overflow; integer arrays have int32 elements. Floats enter as host doubles, must be finite before and after conversion to float, and use binary32 configuration values thereafter. Booleans use host/DS2 boolean binding. Reject wrong types, unknown names and unsupported formats at creation.

| Parameter | Type/default | Phase-1 domain and effect |
| --- | --- | --- |
| clip | Required video node | Fixed-size/format GRAY/YUV/RGB, 1 or 3 planes, UInt8/10/12/14/16 or Float32; positive dimensions and frame count |
| sigma | float 2 | >=0; 8-bit standard deviation; see Wiener specification |
| beta | float 1 | >=1; Wiener gain floor |
| fft_backend | string `pocketfft` | Backend policy in RUN-001; historical default was fftw |
| planes | int array, omitted | Omitted **or empty** selects every actual plane; duplicates have no extra effect; every element in [0,plane_count) |
| bw, bh | int 32,32 | Each >=2; selected-plane geometry must be admitted |
| bt | int 3 | Only 1 implemented; omission still resolves to 3 and fails explicitly |
| ow, oh | int -1,-1 | Any negative means floor(block/3); otherwise 0..floor(block/2) |
| kratio | float 2 | Inactive for bt=1; require default in this batch |
| sharpen | float 0 | Must be 0; nonzero is deferred |
| scutoff, svr | float 0.3,1 | Require defaults while enhancement is deferred |
| smin, smax | float 4,20 | Require defaults |
| measure | bool true | Backend planning hint; no effect for PocketFFT |
| interlaced | bool false | Must be false |
| wintype | int 0 | 0,1,2 |
| pframe, px, py | int 0,0,0 | Require defaults while pattern estimation is deferred |
| pshow | bool false | Must be false |
| pcutoff | float 0.1 | Require default |
| pfactor | float 0 | Must be 0 |
| sigma2, sigma3, sigma4 | float, each inherits sigma | Each must equal parsed sigma; unequal values require deferred frequency-dependent noise |
| degrid | float 1 | >=0; 0 disables, nonzero follows per-block grid subtraction/restoration |
| dehalo | float 0 | Must be 0 |
| hr, ht | float 2,50 | Require defaults |
| l, t, r, b | int 0 each | Must be exactly 0; ROI is deferred |
| opt | int 0 | 0 automatic Highway, 1 own scalar kernels; other values deferred |
| ncpu | int 2 | Positive requested maximum; capability mapping per RUN-001 |
| mt | bool false | Must be false; no internal plane pool in this batch |

Compare default-valued floats after binary32 conversion. A different value on a deferred parameter is rejected even when it would be dormant in the old spatial path. This is an explicit phase boundary, not a claim that the reference rejects it. The legacy negative-ROI clamp and invalid-wintype fallback are not adopted.

## Frame computation

1. At creation resolve format, plane selection and [grid](kernel-geometry.md); build immutable [windows](kernel-window-reconstruction.md) and [noise/grid](kernel-wiener.md) tables for selected planes.
2. Request only clip[n], preserve properties, and copy unselected visible planes.
3. Reflect each selected plane into its cover, subtract its applicable midpoint and gather/window blocks.
4. Forward FFT, Wiener/degrid, normalized inverse FFT and synthesis, using bounded batches without changing contribution order.
5. Crop, quantize/clip once and publish the fully initialized frame.

Execution, stable error ownership, concurrency and frame atomicity follow [RUN-001](../execution.md). Output keeps the source's dimensions, format, count and rate. Finite float samples outside nominal range follow the FFT3D final clipping contract.

## Examples and errors

`core.neo_fft.FFT3D(src, bt=1, fft_backend="pocketfft", opt=1)` exercises the spatial path with historical algorithm defaults, including degrid=1. With src=GRAY8 128x96 the default grid is valid; block size 32 and overlap 10 produce Nx=8,Ny=6, cover 186x142, offset (22,22).

`planes=[]` processes all planes. On YUV420 256x192, `planes=[0]` processes only Y; U/V are byte-for-byte copies. `bt=3`, `sharpen=0.1`, unequal sigma2, nonzero ROI or interlaced=true must fail at creation with an unsupported-feature message. `beta=0`, out-of-range plane index, unsafe tiny-plane padding and overflow are invalid-argument errors. Do not replace a rejected mode with pass-through or with bt=1.

Acceptance must cover [P1-ACCEPT-001](../acceptance.md), including default degrid and float chroma clipping. A valid call is not evidence that its output has been compared to the fixed reference.
