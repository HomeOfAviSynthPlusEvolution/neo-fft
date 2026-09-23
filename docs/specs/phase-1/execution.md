# Phase-1 execution and error contract

Historical phase-1 contract. For current thread controls and removed parameters, apply the later [execution-policy revision](../phase-4/execution.md).

Specification: RUN-001. Evidence: new resource/host contract, preserving stateless spatial behavior.

## Creation and frame dependencies

Creation parses typed parameters, validates finite conversions and active feature gates, resolves plane selection, validates each selected plane's geometry, and constructs immutable configuration/window/grid tables. Check resource sizes before allocation. Initializing a filter does not request video frames. It may fail for invalid parameters, unsupported features, invalid geometry, backend availability or allocation failure.

For output n, request only `clip[n]`. Neither function requests a noise sample, temporal neighbor or another output node in phase 1. Copy source frame properties and preserve dimensions, format, frame count and frame rate. Do not alter color range, chroma location, timing or arbitrary user properties. Unselected planes are bitwise copies of visible source samples. Do not expose internal workspace or FFT metadata as frame properties.

After receiving a frame, revalidate the host descriptor against the plan and acquire exclusive scratch. The source is read-only. Publish the destination only after all selected planes and copies complete. Errors must leave no apparently successful partially written frame.

## Ownership and scheduling

The filter instance owns immutable plans and tables. Each active request owns its accumulator, block buffers, FFT execution resources and temporary spectra. Use RAII to release them on success, failure or cancellation. No process-global numerical state, retained host frame pointers, shared mutable block buffers or RNG is needed.

Phase 1 uses host frame-level parallelism; each request executes synchronously, with one internal worker and one FFT thread. No temporal cache is required. Any bounded workspace pool must permit independent concurrent requests without waiting for a future output request to release a resource. Never hold a pool, cache or FFT planning lock while asking the host for a frame. Do not assume backend execution resources are reentrant merely because the configuration is immutable.

Memory scales as immutable tables plus active requests times workspace, not video length. Bound the number of blocks transformed together by an internal budget. A legal final partial batch must work. Reusing or resizing a workspace cannot retain contributions from a previous plane/frame; zero only the accumulator's full logical extent before reuse.

All repeated requests for n, in one fixed build/dispatch configuration, must return bitwise-identical visible samples regardless of previous frame order and request concurrency. Cross-backend and scalar/Highway comparisons use numerical budgets instead; these are separate claims.

## Execution controls

Both functions use `opt=0` for automatic Highway dispatch and `opt=1` for the project's scalar kernels. This does not claim that the FFT library itself is scalar. Reject other opt values in phase 1 as unsupported execution modes; later compatibility may map legacy ISA selectors explicitly.

Omitted `fft_backend` selects `pocketfft`. Explicit `pocketfft` is required to work. Explicit `fftw` works only if that adapter and runtime are available; otherwise give an error, without fallback. Other strings, including GPU backends, are rejected. The reference defaults to FFTW; differential tests must select the same backend explicitly.

FFT3D `ncpu` remains default 2 as a legacy **requested maximum**. PocketFFT has no internal threading here, so accepted positive requests all resolve to one effective FFT thread; this matches the selected reference backend's capability limitation. `mt=false` is required. `measure` is accepted but has no planning effect for PocketFFT; if FFTW is added, false selects estimate and true selects measured planning on private buffers. Log effective backend/thread/planning choices in the test manifest, not per-frame properties.

DFTTest `threads` and `fft_threads` default to 0 (automatic); nonpositive values mean automatic and resolve to 1. Explicit 1 is supported. Positive values >1 are a phase-1 unsupported execution request, even if a reference accepts or caps them. Any optional FFTW adapter in this batch still executes at one thread; FFT3D explicit ncpu>1 with FFTW must therefore be rejected, and omission resolves to 1 under the documented phase-1 policy.

## Error classes and example

| Stage | Examples | Observable result |
| --- | --- | --- |
| Create | Wrong type, non-finite parameter, out-of-int32 integer, unsupported temporal mode, invalid overlap, padding beyond admitted domain | Filter creation fails with function and offending parameter/domain in message |
| Evaluate | Input descriptor mismatch, selected float plane contains NaN/Inf, non-finite intermediate, allocation/backend failure | Frame request fails; owned error text survives callback return |

Do not compare exact historical error strings or undefined failure ordering. The new messages distinguish invalid arguments from valid-but-unimplemented modes. Validate all supplied floating parameters even if their controlling feature is inactive. Known inactive settings are handled only as documented by each plugin spec.

Example: frame requests [4,0,4,2] require only those source frames and return the same frame-4 samples twice. A failure while processing frame 2 cannot poison frame 4's tables or later accumulator. A copied float plane may retain arbitrary sample bit patterns, while a selected plane containing NaN is a frame error.
