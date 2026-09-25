# Phase-4 acceptance

Specification: P4-ACCEPT-004. These are implementation gates, not claims that tests have already passed. Inherit phase-1/2/3 tests and calibrated error-budget protocol. Keep source/operator correctness, request determinism, ownership and performance claims separate.

## Oracles and evidence

Use independent scalar test oracles for the recurrence, row mapping, diffusion and coordinate hash; do not compute expected values by calling the production helpers. Fix input seeds and preserve the existing calibration seed 1 / holdout 17,101 split for numerical budgets. Freeze budgets by algorithm, format, backend and dispatch before evaluating holdouts. Record commit/compiler/flags, backend/dependency versions, effective worker controls, format/geometry and skipped paths.

Compare old-reference Kalman only in a fresh instance requested sequentially 0,1,..., on finite admitted cases with the phase-3 reference-compatible choices. Kalman sampled comparisons must include multiple positive pfactor values (for example 0.2, 0.7, 1 and 2), which produce identical output for otherwise identical inputs; restricting coverage to pfactor=1 would miss an incorrect inherited Wiener strength multiplier. Reference random-order Kalman and old dither>=2 RNG are not correctness oracles. Zero-noise, unsafe crop/odd-height corners and deliberate phase-3 differences use the independent specification oracle. No unsafe reference execution is needed to prove a rejection case.

## Kalman operator and output

Required direct-kernel cases:

- Zero L and uniform R0 initialization for uniform, analytic and sampled models, including patterned R!=R0. Source frame 0 must not seed the spectrum.
- Both branches; strict equality to motion threshold; only real or only imaginary exceeding it resets both components; kratio=0; first-step examples in the operator; bins at SIMD width-1/width/width+1 and odd transform dimensions.
- Uniform R=0 with X=0 and X!=0 stays finite; pattern P=0 uses 1e-15; all positive sampled pfactor values select the same unscaled power before the floor, including finite strengths whose irrelevant multiplication would overflow or underflow; no temporal multiplier, Wiener beta or pre-recurrence degrid.
- Active overflow fails cleanly; large squared differences that trigger reset do not fail merely because the comparison temporary is infinite. Masked SIMD lanes do not leak speculative NaNs. State/new-state aliasing follows the production kernel's documented contract and guard buffers show no tail overwrite.
- Finite large sigma that is inactive for phase-3 sampled Wiener becomes active initial covariance in Kalman and is rejected if R0 overflows; raw invalid kratio fails even in other modes, but inactive kratio-derived overflow does not.
- Sharpen/dehalo affect output only. Compare stored L,C,Q with enhancement enabled/disabled over several steps. Intermediate replay must not feed clipped or enhanced results back.

Required full-filter cases include all supported sample formats and GRAY/YUV/RGB, selected/subsampled planes, normal/analytic/sampled models, ROI and interlaced combinations. Frame 0 and a one-frame clip are bitwise copies, including non-finite float pixels, with unchanged properties; n>0 rejects only consumed non-finite source regions. Test sampled pframe=0 and pframe beyond target n, and all inherited pframe clamping rules.

Recurrence can amplify a tiny numerical change into a different motion branch. Validate direct state transitions against their oracle, and freeze separate multi-step scalar/Highway output/state budgets on long constant, smooth-motion, threshold-adjacent and scene-cut sequences. Do not use a mean error alone to conceal wrong branch/reset behavior. Same-arithmetic execution-order comparisons are always bitwise.

## Replay, staged release and cache

Revised 2026-09-24: compare bounded cold output against an independent recurrence over max(1,n-W)..n. Compare checkpoint continuation against the leased state's chronological descendants, not a universal canonical output. Cover W=0/4/8/16/INT32_MAX, near-start clamp, n=10000, int32 end, earliest eligible and just-too-old checkpoints, exact/nearest hits, backward seeks, disabled cache and eviction. Retain a snapshot while evicting it; verify first retained publication wins at duplicate n without shared mutation. Concurrent same-history descendants agree; different histories need not. Full-prefix canonical tests must explicitly configure enough warmup or stay within the window.

Instrument the staged provider/VS bridge, not just core scratch allocations. Verify at most one acquired source-frame owner retained at replay stage boundaries per request and bounded Neo-FFT/DS2-owned storage for fixed geometry/concurrency. Independently count VS frameCtx-held pixel references; fake-host tests must model this second ownership layer and aliased input nodes. Test the actual target runtime over increasing replay distances before claiming bounded dependency pixels; ordinary staged integration tests and API presence are insufficient. Record host/version and report host history metadata, upstream caches and property-held references separately, without an unconditional total-host-memory bound. Check release then lookup error, release then reacquire same identity, pending-frame release rejection, pframe overlap/reacquisition, final n property copy, and the unchanged behavior of consumers that never call release. Owning get() copies must remain valid across store release and disappear when their last owner expires. Confirm all request-owned references are released after completion, cancellation and injected stage failure.

Use a tiny-plane synthetic provider for large logical frame numbers. Verify cold dependencies remain within W+1 recurrence frames, even near INT32_MAX, and that W=INT32_MAX never overflows arithmetic. Exercise long staged-release paths by explicitly increasing W; default cold work must not grow with n. Model pframe acquisition is separate from the recurrence budget and requires its own source-closure check.

Inject failures at sample acquisition/build, intermediate source delivery/FFT, checkpoint clone, final allocation/conversion and worker execution. Failed requests publish no partial checkpoint; a retry and an unrelated concurrent request match clean results when initialized from the same state and consuming the same history. Ensure no host fetch/wait occurs under cache/model/pool locks and no request depends on an unfinished output of the same filter. Capture the actual new pinned DS2 revision and passing DS2 lifetime/integration tests.

## ROI and field packing

- Direct row labels prove H=2,4,8 mappings and exact inverse. Distinguish full-height packed filtering from two independent field filters with a signal crossing the middle seam.
- Nonzero unequal l,t,r,b, ROI-relative odd top origins, aligned 420/422 chroma, RGB and luma-only selection. Reject selected-chroma misalignment, negative/empty crop, unsafe reflection, and odd selected packed height. Unselected chroma does not add an alignment/even-height restriction.
- Verify copied exterior/unselected pixels byte-for-byte, including float NaNs outside consumed regions; metadata always comes from n. Guard pitches and row endpoints.
- Exercise bt=-1,0,1,2,3,4,5; first/last temporal fallback; sampled pattern/manual/automatic selection and pshow in cropped/packed coordinates. No temporal index doubling or property-based field-order switch.
- Differential tests on safely aligned even-height reference inputs. For nonaligned crops/odd packed heights, assert this specification's rejection without running reference UB.

## Dither

Direct oracle inputs: width 1,2 and SIMD-adjacent widths, height 1/2, fractional ramps, negative/out-of-range finite values and saturation near 0/255. Distinguish historical E-D from corrected-error alternatives, clipped D from unclipped quantization and discarded edge weights from renormalization. Use values where error crosses multiple row/block boundaries; no strip resets or serpentine scans.

Check all hash golden vectors exactly with an independent unsigned implementation. Exercise omitted seed versus 0, seeds 1/17/INT_MAX, dither 0/1/2/large positive/INT_MAX and negative rejection. Confirm mode 1 ignores seed, mode 0 matches inherited conversion, and modes >=2 vary with coordinates/seed while remaining reproducible. An integer-overflow-free scale computation is mandatory.

Integration covers both smode values, temporal T=1 and odd T>1, all five ftypes, zmean on/off, curves and sampled nlocation. Dither is applied once after final accumulation/crop. Explicit planes=[] and unselected planes copy exactly; high-depth/float dither>0 equals their dither=0 output. Reordering/duplicating planes does not alter noise identity. Mode 1 has a reference comparison; randomized modes have the specified-hash oracle.

## Execution and release gate

Validate all consumed public control domains and historical opt aliases, including kalman_warmup in inactive modes. Verify signature-only mt/ncpu/measure/threads/fft_threads/fft_backend values are not read or semantically validated by the plugin; registered host types still apply. Both filters execute on the calling host worker, with PocketFFT internally single-threaded. Ignored values cannot alter workspaces or output. KernelInfo retains the fft_backend diagnostic but omits threads and fft_threads. Stateless modes and Kalman requests with identical initial state and consumed history retain bitwise results under concurrent host requests. Scalar versus Highway uses calibrated budgets. Confirm one instance cannot alter another instance's dispatch.

Stress active concurrency within the test machine's configured limit (build/test processes together <=4), repeated create/destroy, warm/cold model/cache, error retries and cancellation. Check idle workspace/checkpoint budgets separately from active leases and disclose untracked backend allocations. Use supported sanitizer/guard/ownership instrumentation; list unsupported or skipped checks rather than treating their absence as a pass. Error assertions use structural stage/category, never localized text matching.

Final evidence must contain: DS2 revision `f1d51bd0217d3878f995375e95c2827b4604facf` actually selected by the build, independent operator results, admitted-reference comparisons, VS end-to-end results, bounded Kalman request-history matrix and deterministic stateless-mode matrix, bounded-retention counters, failure/lifecycle results, inherited regressions and declared optional skips. Check that an existing CMake cache or FetchContent source override has not silently selected the old dependency. A benchmark or successful compilation alone does not satisfy this gate. No production AVS/GPU/FFTW claim follows from the required VS PocketFFT acceptance.
