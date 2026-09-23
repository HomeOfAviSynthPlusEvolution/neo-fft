# Phase-4 CPU execution and lifetime

Specification: RUN-004. Extends [phase-3 execution](../phase-3/execution.md). Kalman dependencies use [replay](replay.md); all other modes keep their declared closures. This document opens execution controls without changing mathematical modes.

## Public control mapping

| Function/control | Parsing and effective behavior |
| --- | --- |
| FFT3D opt | Any int32 accepted, matching the pinned dispatcher; 1 selects own scalar kernels, every other value automatic Highway |
| DFTTest opt | Only 0,1,2,3,8; 1 selects own scalar, all others automatic Highway |
| FFT3D mt | Removed by the user-directed execution-policy revision; supplying it fails creation. All plane work runs on the calling host worker |
| FFT3D cache_frames | Default -1: auto `bt + host_threads - 1`; 0 disables raw-frequency caching; positive int32 is a frame limit; values below -1 fail |
| FFT3D cache_mb | Default 128 MiB; -1 selects the same default budget; 0 disables raw-frequency caching; positive int32 is a memory limit; values below -1 fail |
| DFTTest threads | Reserved, default 0; any int32 accepted, <=0 normalizes to 1 and positive values clamp to 16; currently has no effect on execution or workspace retention |
| fft_backend | PocketFFT required/default; unavailable explicit FFTW fails creation; no silent fallback |

FFT3D ncpu and measure, and DFTTest fft_threads, have been removed by the user-directed interface cleanup. Supplying them fails creation, including historical defaults. No requested FFT-worker state is retained. KernelInfo's read-only fft_threads=1 still reports the backend's effective count; it is not a filter parameter.

opt aliases do not force a named instruction set or bypass runtime CPU capability checks. opt=1 controls own kernels, not the third-party FFT's ISA. Do not accidentally map historical opt=8 to a GPU backend.

FFT3D has no internal executor, including ordinary filtering, preview, Kalman replay and rendering. Process planes in order on the calling host worker. Host frame concurrency and cache/model synchronization remain supported. This user-directed revision supersedes the historical mt parameter and plane-worker requirements in earlier phases.

DFTTest has no internal executor. Block transforms/filtering execute sequentially in bounded batches on the calling host worker. The threads parameter retains its registration, default and parsing rules as a future reservation; it does not create workers or size workspace retention. Host requests may still run concurrently with private workspace leases. This user-directed revision supersedes earlier own-worker requirements, including phase 5. Do not block a host callback waiting for a workspace owned by an unrelated request that itself needs host progress.

PocketFFT remains single-threaded internally and has no user-controlled measured planning. This revision supersedes earlier contracts for ncpu, measure and fft_threads, including the optional FFTW control mapping; it does not implement another backend. Record the effective execution settings in test manifests; no extra frame properties are added.

## Deterministic work partitioning

Window/model construction order and scalar operation definitions remain inherited. Sequential block processing preserves fixed FFT groups and commits overlap-add contributions in the same canonical order as the scalar path; racing float atomics, independent strip sums merged in a different order or scheduler-dependent reductions are forbidden. Keep batches bounded independently of image frame count and replay distance. A scheduling change cannot alter padding, logical bins, tail safety or reflection.

Kalman time steps and selected planes execute sequentially within each request; independent host requests retain private mutable state. Dither's final scan runs sequentially per plane. Pattern/noise estimates use a fixed reduction order and immutable publication.

Within one build/backend/dispatch/plan arithmetic configuration, output must be bit-identical across request order, host concurrency, reserved threads values, cache warmth, replay partition and repeated runs. Disabling/enabling checkpoint caching cannot change it. Scalar versus Highway and distinct FFT plans/backends use the inherited calibrated numerical budgets; do not silently require bit equality between different FFTW threading plans. Dither random samples themselves are exact and independent of every execution control. Test near quantization/motion thresholds, where small spectral differences can change discrete decisions; report these separately instead of averaging them away.

## Resource accounting and errors

Immutable windows, models, plans and checkpoint shape belong to the instance; host frames, removed means, spectra, inverse/accumulation buffers, dither rows and mutable Kalman state belong to a request/leased workspace. No scratch array is shared concurrently merely because two requests have the same frame index or worker id. Worker id is not RNG identity.

Count sizes, pitch/row offsets, grid products, FFT strides, state arrays, model arrays, batch scratch and worker multiplicities with checked arithmetic before allocation/narrowing. Large valid requests may fail cleanly for allocation; they cannot wrap into smaller buffers. No work queue/cache/list has capacity proportional to clip length. Dither error rows require 2*W*sizeof(float) per active UInt8 plane, not a full error history.

Idle workspace retention is capped at 64 MiB and one reusable workspace object per instance across all planes, independently of the reserved threads value; count actual capacities including pooled dither/FFT buffers. A larger returned workspace is freed instead of cached. Active requests may allocate outside the idle pool and must be included in peak-memory diagnostics; the idle budget is not a global process-memory guarantee. The separate checkpoint budget is defined in replay.md. Do not retain request-bound host-frame owners in an idle workspace. Backend/library allocations must be disclosed separately if they cannot be measured by the pool counter.

No exception crosses a C/host callback boundary. Preserve an owned error string for as long as the host/DS2 can use it; never return a pointer into a destroyed exception or temporary string. Identify invalid-argument, unsupported-feature, resource and frame-processing failures structurally in tests, by stage/category; do not infer outcomes from localized message substrings. Human messages should include the function, offending parameter or frame/plane and meaningful context without promising identical reference wording.

Only publish fully initialized output. On any processing failure discard that output/private state and return one stable error; do not return a partly copied frame. Other in-flight requests retain valid model/checkpoint/workspace leases. Release handles exactly once on success, error and cancellation. A failed speculative duplicate model build cannot invalidate a ready published model.

## Optional FFT3D input-spectrum cache

For ordinary bt=2..5 requests, an instance may share complete finite raw spatial FFT spectra, after input centering, reflection and analysis windows but before any grid subtraction, filtering or enhancement. A key is the actual source frame, selected plane and spatial FFT block row within one immutable source/configuration/ROI/backend instance. Packed interlaced geometry is part of that instance. Endpoint single-frame fallback may reuse the same raw product. bt=-1/0/1 and preview bypass this cache.

Both public limits apply together, allocation is demand-driven, and no limit depends on clip length. Auto frames uses AviSynth's CACHE_INFORM_NUM_THREADS when delivered, otherwise one worker; VapourSynth resolves from core.numThreads at creation. This is a capacity heuristic for adjacent requests, not a guarantee about actual in-flight frames. Positive explicit frame limits are not changed by thread notifications.

Each retained block row reserves only its own spectrum payload plus actual entry/control/map-node storage. Distinct retained source frames count toward the frame limit and their map-node storage is also charged. Requests lease only their current block row across temporal slots and release those leases before advancing to the next row. For registered requests, unpinned rows with no remaining registered consumer are evicted first; ties prefer earlier source frames, then plane/row order. Frame-count pressure similarly prefers a fully unpinned frame with no remaining consumer of its retained rows, then earlier frames. Remaining future consumers are only a soft priority: when all candidates are needed, eviction still proceeds within the hard budget. Callers that never register retain the former LRU fallback. The memory limit includes building and leased entries; pinned entries cannot be evicted or uncharged. No spectrum payload is allocated before admission. A small candidate bookkeeping allocation may precede admission; the fixed cache object, allocator overhead and request workspaces remain outside this entry budget. A miss that cannot reserve capacity uses the existing bounded-block path without waiting for memory. Decreasing an auto frame limit trims unpinned entries; existing pins remain valid and are trimmed on subsequent access after release.

Ordinary host requests register their actual source-frame interval before requesting input dependencies, including the single-frame endpoint fallback. Each plane records the first unfinished block row; completing a row removes that request's demand for earlier rows. Registration ends with processing, exception or host cancellation through request-state RAII. Concurrent duplicate requests have separate registrations, and plane progress may advance independently. The cache must outlive its registrations; host request teardown precedes instance destruction. Registration bookkeeping belongs to active request state, scales with outstanding requests rather than clip length, and is outside the retained-entry budget. It never owns a host frame or spectrum lease. Earlier-frame eviction favors forward playback; undeclared future or reverse requests may recompute evicted spectra. Registration is an eviction hint, never a correctness prerequisite or a reason to wait for memory.

Same-frame/plane/row construction is single-flight. Waiters release the cache mutex. The producer executes synchronously from already acquired source pixels, does not call the host, acquire another cache entry, or wait for another output request, and either publishes a complete finite row or signals failure. Thus there is no request dependency cycle. Failed products are never returned; existing waiters receive the failure, future requests may recompute uncached until eviction. Successfully built independent raw rows can survive a later output failure. They contain no mutable filter output and no host frame owners.

The cache must preserve bit-identical output within one arithmetic configuration for enabled/disabled, small-budget fallback, eviction, sequential/reverse/shuffled/repeated and concurrent requests. Memory accounting and publication/failure tests are required independently of timing. Instance-local diagnostic counters expose row requests, reuse hits (including joined builders), waits, builds, capacity/failure bypasses, evictions and requested/computed bins, active/peak request registrations, and evicted rows/consumer counts that still had registered future demand. Miss rate counts builds plus bypasses; compulsory new-input transforms are misses too. These counters do not add host properties or a public plugin API.

## Teardown and dependency gate

The implementation must make instance ownership outlive its active host requests and workspace leases before destroying plans, models and allocation providers. There are no internally owned workers to join. Never call a host environment after its lifetime ends. No cache lock is held while releasing a host frame through a callback that could reenter.

The DS2 release facility in replay.md is mandatory before enabling production Kalman replay. Integrate the committed/pushed revision `f1d51bd0217d3878f995375e95c2827b4604facf` through `cmake/Dependencies.cmake`, retain feature/configuration checks and run staged ownership tests on the actual linked version. No additional DS2 release-interface implementation is required. Separately verify that the target VS runtime effectively releases frameCtx-held dependency pixels; record supported runtime conditions and account for host history metadata/caches separately. Do not claim bounded source memory based solely on one-frame stages, store-owner counts or a non-null releaseFrameEarly pointer.

Phase-4 completion requires the supported VS CPU modes to coexist in one build with inherited phase-1/2/3 acceptance. It does not require AVS entry points, optional FFTW, GPU or unsupported tmode=1. Record skipped optional paths explicitly.
