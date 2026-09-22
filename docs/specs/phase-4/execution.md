# Phase-4 CPU execution and lifetime

Specification: RUN-004. Extends [phase-3 execution](../phase-3/execution.md). Kalman dependencies use [replay](replay.md); all other modes keep their declared closures. This document opens execution controls without changing mathematical modes.

## Public control mapping

| Function/control | Parsing and effective behavior |
| --- | --- |
| FFT3D opt | Any int32 accepted, matching the pinned dispatcher; 1 selects own scalar kernels, every other value automatic Highway |
| DFTTest opt | Only 0,1,2,3,8; 1 selects own scalar, all others automatic Highway |
| FFT3D mt | false: own work runs on the calling host worker; true: selected planes may run concurrently, maximum actual plane count (<=3) |
| FFT3D ncpu | Positive int32 requested maximum FFT workers; default 2; PocketFFT effective 1 |
| DFTTest threads | <=0 resolves to 1; positive values clamp to 16; requested maximum own workers, further limited by available independent work |
| DFTTest fft_threads | <=0 resolves to 1; positive int32 requested maximum FFT workers; PocketFFT effective 1 |
| FFT3D measure | Planning hint only; PocketFFT ignores it; never changes semantic mode |
| fft_backend | PocketFFT required/default; unavailable explicit FFTW fails creation; no silent fallback |

opt aliases do not force a named instruction set or bypass runtime CPU capability checks. opt=1 controls own kernels, not the third-party FFT's ISA. Do not accidentally map historical opt=8 to a GPU backend.

Implement the opened own-worker paths: FFT3D mt=true can schedule independent selected planes; DFTTest threads>1 can schedule independent block transforms/filtering within bounded batches. One worker is valid for insufficient work or resource contention, but accepting the parameter while unconditionally ignoring it is not completion. Test with enough independent work and an internal active-worker counter. No throughput target is imposed by acceptance.

Automatic own-worker policy is fixed at one for DFTTest; FFT3D mt=false likewise adds no own workers. The user may opt into internal work alongside host-level frame concurrency. Share each instance's worker capacity across its requests rather than creating a new persistent pool per frame. Cap additional active own workers to the resolved maximum (FFT3D <=3, DFTTest <=16); the calling host thread may participate. A caller that cannot obtain an internal worker can perform work inline. Do not block a host callback merely waiting for an idle workspace/worker owned by an unrelated request that itself needs host progress.

Required PocketFFT remains single-threaded internally. If optional FFTW is delivered, report a checked effective FFT-worker maximum bounded by 16 and available backend capability; unsupported threaded capability resolves to 1. To avoid nested internal teams, choose own parallelism or FFT parallelism for a request: if effective FFT workers>1, effective own workers=1. Plan creation/destruction and any global backend thread/planner state must be synchronized without affecting active plans. Changing one instance's controls must not mutate another's plan. Record requested/effective controls in test manifests; no extra frame properties are added.

## Deterministic work partitioning

Window/model construction order and scalar operation definitions remain inherited. Parallel block processing stores private results and commits overlap-add contributions in the same canonical order as the scalar path; racing float atomics, independent strip sums merged in a different order or scheduler-dependent reductions are forbidden. Keep batches bounded independently of image frame count and replay distance. A scheduling change cannot alter padding, logical bins, tail safety or reflection.

Kalman bins/planes within one time step may run concurrently; the next source time step begins only after all selected-plane state updates complete. Each request has private mutable state. Dither's final scan runs sequentially per plane; planes are independent. Pattern/noise estimates use a fixed reduction order and immutable publication.

Within one build/backend/dispatch/plan arithmetic configuration, output must be bit-identical across request order, host concurrency, own-worker count, cache warmth, replay partition and repeated runs. Disabling/enabling checkpoint caching cannot change it. Scalar versus Highway and distinct FFT plans/backends use the inherited calibrated numerical budgets; do not silently require bit equality between different FFTW threading plans. Dither random samples themselves are exact and independent of every execution control. Test near quantization/motion thresholds, where small spectral differences can change discrete decisions; report these separately instead of averaging them away.

## Resource accounting and errors

Immutable windows, models, plans and checkpoint shape belong to the instance; host frames, removed means, spectra, inverse/accumulation buffers, dither rows and mutable Kalman state belong to a request/leased workspace. No scratch array is shared concurrently merely because two requests have the same frame index or worker id. Worker id is not RNG identity.

Count sizes, pitch/row offsets, grid products, FFT strides, state arrays, model arrays, batch scratch and worker multiplicities with checked arithmetic before allocation/narrowing. Large valid requests may fail cleanly for allocation; they cannot wrap into smaller buffers. No work queue/cache/list has capacity proportional to clip length. Dither error rows require 2*W*sizeof(float) per active UInt8 plane, not a full error history.

Idle workspace retention is capped at 64 MiB and at most the resolved own-worker count of reusable workspace objects per instance; count actual capacities including pooled dither/FFT buffers. A larger returned workspace is freed instead of cached. Active requests may allocate outside the idle pool and must be included in peak-memory diagnostics; the idle budget is not a global process-memory guarantee. The separate checkpoint budget is defined in replay.md. Do not retain request-bound host-frame owners in an idle workspace. Backend/library allocations must be disclosed separately if they cannot be measured by the pool counter.

No exception crosses a C/host callback boundary. Preserve an owned error string for as long as the host/DS2 can use it; never return a pointer into a destroyed exception or temporary string. Identify invalid-argument, unsupported-feature, resource and frame-processing failures structurally in tests, by stage/category; do not infer outcomes from localized message substrings. Human messages should include the function, offending parameter or frame/plane and meaningful context without promising identical reference wording.

Only publish fully initialized output. On any worker failure join/cancel outstanding work safely, discard that output/private state and return one stable error; do not return a partly copied frame. Other in-flight requests retain valid model/checkpoint/workspace leases. Release handles exactly once on success, error and cancellation. A failed speculative duplicate model build cannot invalidate a ready published model.

## Teardown and dependency gate

The implementation must make instance ownership outlive its active requests and workers, then stop/join its worker resources before destroying plans, models and allocation providers. Never call a host environment after its lifetime ends. No cache lock is held while joining a worker or releasing a host frame through a callback that could reenter.

The DS2 release facility in replay.md is mandatory before enabling production Kalman replay. Integrate the committed/pushed revision `f1d51bd0217d3878f995375e95c2827b4604facf` through `cmake/Dependencies.cmake`, retain feature/configuration checks and run staged ownership tests on the actual linked version. No additional DS2 release-interface implementation is required. Separately verify that the target VS runtime effectively releases frameCtx-held dependency pixels; record supported runtime conditions and account for host history metadata/caches separately. Do not claim bounded source memory based solely on one-frame stages, store-owner counts or a non-null releaseFrameEarly pointer.

Phase-4 completion requires the supported VS CPU modes to coexist in one build with inherited phase-1/2/3 acceptance. It does not require AVS entry points, optional FFTW, GPU or unsupported tmode=1. Record skipped optional paths explicitly.
