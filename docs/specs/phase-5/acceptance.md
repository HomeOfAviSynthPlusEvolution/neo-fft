# Phase-5 acceptance

Specification: P5-ACCEPT-005. Required: independent geometry/operator oracles, full new VS path versus the fixed original AVS reference, own scalar/Highway verification, and regression of phases 1-4. This specification does not claim those runtime gates have passed.

## Geometry and window oracles

Independently implement both the closed-form start range and the original cold GetFrame_T start loops, retaining only blocks containing n. Exhaust all admitted T=1..15/O pairs over first/last frames and several interior strides, odd/even clip lengths, N=T and N>T. Compare starts, target slices, clamp maps and the <=2*T-1 unique-source bound. Include all table examples, negative starts, O=0, exactly half overlap and valid/invalid heavy overlap. A cold ring may compute extra earlier blocks; compare contributing blocks only.

Near integer limits use the direct formulas and a local enumeration of nearby lattice points, not a billions-step original loop or enormous video. Prove floor-division behavior for negative numerators and checked host narrowing. A far-index request's geometry work must stay bounded by ceil(T/(T-O)).

For all twelve temporal windows, odd/even T and all admitted O, verify raw samples, denominator accumulation order, normalized samples, stored h and wscale. Test non-finite/zero phase-energy rejection, zmean's DC admission, zero-overlap sign preservation and singleton behavior. Identity-oracle gain uses stored float h. Freeze inherited calibration seed 1 / holdout 17,101 numerical budgets before holdout comparisons.

## Operator and parameter matrix

- Both spatial modes; odd/even spatial transform sizes where admitted; non-square frames, spatial padding boundaries and subsampled planes.
- All five ftypes, scalar/curved/axis profiles, active nlocation and overridden unused arrays, zmean on/off, negative float chroma, default and nondefault window parameters.
- Full-volume zmean versus mistaken per-slice means; the small temporal-mean and irregular-overlap examples in the operator. Even-T temporal Nyquist bins and SIMD tails with guard buffers.
- sigma/PSD energy recalibration after temporal normalization; raw h2 sampling window and wscale2/wscale ratio with O changed. Sample start not on the output lattice; duplicate sample tuples; first/last valid sample intervals and invalid fn+T.
- smode=0 adds all temporal contributions rather than overwriting; smode=1 preserves time/start/spatial addition order. No per-block rounding or final weight correction.
- dither 0,1,>=2 and seeded reproducibility after full temporal accumulation. UInt8 versus high-depth/float activation remains phase 4. Plane order/duplicates, explicit planes=[] and copied unselected planes/properties.
- All invalid mode/T/O/heavy-overlap combinations, raw invalid values even when inactive, and failures only for consumed derived arithmetic. tmode=0 continues ignoring supplied tosize after int32 parsing, exactly as previously specified.

## Original-reference differential

Build original pinterf/dfttest at `465d5d184b1be47244b7f8b4f3da7a0d57aa4e29`; record compiler, FFTW binary/version, plugin hash and original source revision. Use its scalar opt=1, threads=1 and dither=0 as the baseline. Generate identical visible samples in AVS and VS and compare exported planes/frame indices, not container encoding or host-dependent properties. Explicitly map all differing defaults: new sigma=8, tbsize=3, sbsize=16,sosize=12, and any other varied control. Disable original lsb/lsb_in modes. Use planar common formats; map planes to the original Y/U/V flags by its actual registration semantics.

Mandatory reference cases cover tmode=1 with T=2/3/4/5/6, O=0, partial, half and heavy overlap, both spatial modes, all ftypes, all admitted temporal windows, zmean, endpoints and interior block-phase transitions. Include UInt8, supported high-depth integer and float, plus subsampled YUV and planar RGB. If a format or combination cannot run in the chosen original binary, record it as uncovered and use the independent oracle; do not silently mark full original-format equivalence. Before final acceptance resolve missing mandatory common-format cases.

Adapt numeric curve arrays/noise tuples to the original string parameters using enough digits to preserve binary32 values. Compare only semantically matching configurations. Retain phase-3 declared curve/singleton/sampling differences; isolate them from new temporal behavior. Phase-4 randomized dither is compared with its own hash oracle, never with the original evolving RNG. Dither=1 may receive a separate matched-conversion comparison; core temporal differential uses dither=0 to expose floating errors before diffusion.

Request original outputs in sequential and fresh-instance random order to characterize its ring reuse. The independent mathematical contract remains authoritative if an original cache bug/order-dependent result is found: isolate the case and describe the deviation instead of importing mutable-state behavior. Existing neo-dfttest is not a tmode=1 oracle because it rejects that mode. Original binary execution being unavailable/skipped leaves this differential gate open; scalar-only agreement is not a replacement.

New VS properties are checked against source n, including when n is not the first requested index. Native host-specific property differences are not algorithm pixel errors. Phase 6 remains responsible for Neo-FFT's own AVS entry points and two-entry equivalence.

## Determinism, resources and regression

Within one build/backend/dispatch/arithmetic configuration compare bitwise: sequential, reverse, shuffled, repeated, concurrent duplicate and disjoint targets; own-worker counts; cold/warm models; optional cache enabled/disabled, eviction and duplicate publication. Verify canonical contribution order at overlap/quantization thresholds. Scalar versus Highway/different FFT backends use separately frozen budgets and disclose branch-sensitive differences.

Inject source/FFT/allocation/model/output errors and cancellation. Verify retry matches clean results, no partial frame/model/cache publication, no leaked owner/lease and no shared mutable output accumulation. Instrument declared dependency union, actual live source owners, bounded batch/cache/workspace bytes and lifecycle cleanup. Vary clip length/target index at fixed geometry without allocating a clip-length fixture; state must not grow with index. Benchmark separately from correctness.

Rerun phase-1/2/3/4 accepted subsets. Supersede their old tmode=1 rejection assertions with this admitted-domain matrix; preserve all remaining historical tests, especially tmode=0 odd parity, its raw-window calibration and boundary slots, FFT3D Kalman/ROI and deterministic dither. No FFT3D numerical change is expected.

Build/test total concurrency remains <=4 and configurations run sequentially. Record supported sanitizer/guard checks and optional backend/platform skips. Closure requires the new VS public path, independent oracles, pinned original differential, Highway comparison, deterministic scheduling/resource evidence and inherited regressions; successful compile or a benchmark is insufficient.
