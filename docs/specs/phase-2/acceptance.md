# Phase-2 acceptance

Specification: P2-ACCEPT-001. This is a protocol, not a passing result. Extend [P1-ACCEPT-001](../phase-1/acceptance.md), retaining pinned binary provenance, public VS invocation, deterministic fixture generation, calibration/holdout separation, error categories and per-plane metrics. Phase 2 stays open if required runs are unavailable or skipped.

## Comparison and numerical budgets

Run the phase-1 A/B/C/D pairs with the temporal configurations below: pinned reference scalar versus new scalar, pinned scalar versus new Highway, pinned Highway versus new Highway, and new scalar versus new Highway. Record actual reference registrations, DLL hashes, commits/dirty state, build recipe, dependencies, host/Python, backend, ISA and effective thread settings. Old released DLLs are supplemental, not replacements for the pinned modernized source. An identity or new-scalar/new-Highway-only test cannot establish reference equivalence.

Provisional outer targets are float32 maximum absolute error <=2e-5 in native float sample units and integer maximum difference <=1 LSB; these are **not experimentally established backend budgets**. Report maximum absolute per-sample error, MAE, RMSE, differing fraction and coordinates separately. Before holdout tests, freeze finite per-algorithm/format/backend/dispatch budgets under phase-1 rules; any chosen limits must also meet the provisional outer targets unless a documented numerical review revises them. No missing budget, NaN, systematic scale error or wrong branch passes by tolerance. Never loosen limits merely to accept current output.

Use all phase-1 formats/plane/property baseline cases with bt=3 or T=3 and clips N>=3; run detailed interactions on Gray float and high-bit YUV, plus pairwise combinations. Use zeros, constants, distinct per-frame constants, temporal ramps, impulses in each temporal slot, alternating frames, spatial sinusoids/checkers and seeded noise. Supply N>=15 cases. Keep deterministic seeds 1 for calibration and 17/101 for holdout.

## Mandatory matrix

| Case ID | Required cases | Checks |
| --- | --- | --- |
| P2-API | Omitted bt/tbsize; explicit 1; FFT3D 1..5; DFTTest 1,3,5,7,9,11,13,15 | Default 3 now works on sufficiently long clips; parity is algorithm-specific |
| P2-REJECT | DFTTest T=0,2,4,14,16,17,T>N; tmode=1; nonempty noise/profile arrays; non-finite/overflow values; FFT3D bt=-1,0,6 | Correct creation stage/category; no silent coercion, including DFTTest planes=[] |
| P2-NORMALIZE | DFTTest tmode=0,tosize=0,-1,INT32_MAX; out-of-int32 argument | In-range values become effective 0 before overlap validation; narrowing overflow rejected |
| P2-F3D-TIME | bt=2..5, sigma=0,2,8, beta=1,2, degrid=0,0.5,1 | Correct cur phase, nontrivial attenuation, PSD rather than magnitude, T-scaled noise |
| P2-F3D-GRID | Different DC per frame and spatial block; impulse in one neighbor; nonrectangular spatial windows | Grid from cur only; T*M removed from temporal DC; no per-neighbor degrid |
| P2-F3D-EDGE | Every frame for N=1..7 and every bt=2..5 | Exact fallback predicate and explicit bt=1 comparison; bt=2 last frame remains temporal where legal |
| P2-DFT-SLOTS | Every edge slot for T=3,5,15; N=T and N>T; inspect dependency log | Exact slot count/order and endpoint replication; no compaction or extra requests |
| P2-DFT-FILTER | Every ftype=0..4, zmean true/false, smode 0/1, T=3,5; nontrivial phase-1 filter values | Power/threshold units, sigma2 branch, mean restoration; center synthesis |
| P2-DFT-BETA | T=3, f0beta=0.5,1,2 and float neighbors at dispatch boundaries | Strict exponent dispatch thresholds, not only default ftype=0 |
| P2-DFT-WIN | All twin=0..11 at T=1,3,5,15; Kaiser beta 0,2.5; spatial windows 0,6,7,9 and both spatial modes | Stored-float E/wscale, actual non-unit center gain, asymmetric window 9, no temporal normalization |
| P2-DFT-MEAN | Constant volume; q=[1,2,6] example; zmean on/off; zero-template-DC case | One full-volume ratio, residual before PSD, restoration before inverse; controlled zero-DC creation error |
| P2-DFT-GEOMETRY | smode=1 B=8,9,16 with zero/half/legal heavy overlap; smode=0 B=1,3,9 | Correct common spatial origins across time, odd/even compressed width, crop and tails |
| P2-COPY-META | All/omitted/empty/subset/duplicate planes; per-frame sentinel properties | Original differing empty-list behavior; bitwise copies and properties from n |
| P2-REGRESSION | Phase-1 supported bt=1/T=1 cases and deferred-mode rejections | Earlier numerical behavior preserved; replace only old 'temporal size unsupported' assertions |

Invalid window/mean combinations belong in rejection coverage, not successful numerical coverage. Ordinary PSD filtering is in scope; **nlocation noise estimation is not**. Verify nonempty nlocation is rejected; do not require profile sampling or silently implement a subset in this phase.

## Independent numerical and storage checks

1. Implement binary64 direct DFT checks for [FFT-TIME-001](kernel-fft.md): forward sign/axis order, full Hermitian pairing, Parseval, 1/V inverse, non-square shapes, strides and active batches. Round-trip alone misses matching forward/inverse errors.
2. For FFT3D compare chronological inverse-at-c with current-at-zero unrolled evaluation at T=2..5. Use differing complex slots and nonzero noise/degrid as well as identity. Assert current-only grid ratio and noise=((T*sigma_eff)*sigma_eff)/norm.
3. For DFTTest independently build float h, sum E in z/y/x order and compare wscale and A/B2/L/H, including the rectangular/Hann examples in [DFT-TIME-001](dfttest/kernel-temporal.md). Compare spectrum/template/residual powers and pre-quantization center contributions; pixel-only checks can conceal a scale error.
4. Inject R=3+4i (p=25) into the **filter kernel**, with exact/adjacent thresholds, both type-3 endpoints, zero residual and type-4 H=0. Test f0beta boundaries. Arbitrary complex kernel fixtures need not be valid C2R volumes; complete inverse fixtures must preserve Hermitian pairing.
5. Check zmean with full-volume constants and changing temporal constants. For B=1,T=3, rectangular windows, ftype=2,sigma=0,zmean=true, q=[1,2,6] reconstructs 3 in 8-bit amplitude units. A per-frame mean implementation returns 2 and fails.
6. Check synthesis uses (U*V)*h at center, full-volume energy, no forced unity and no extra 1/T. Validate tw=6 non-unit gain against a double analytic/table oracle before output quantization.
7. Guard natural-but-not-vector-aligned rows, row/volume/batch padding, active counts 0,1,C-1,C, odd widths and SIMD tails; test checked dimension/stride/byte overflow and failed allocation cleanup. Selected neighbor NaN/Inf must fail; copied unselected samples retain bit patterns.

## Scheduling, concurrency and closure

Use an instrumented provider or host fixture to record declared/requested source frame indices and slot reuse. Verify all DFTTest duplicates remain numerical slots, fallback makes no unused-neighbor requests, and properties come from n. Inject upstream failures and confirm a later request can reuse a workspace safely.

For a fixed build/backend/dispatch compare sequential, reverse, seeded shuffled and repeated requests, simultaneous duplicate n, overlapping neighborhoods and separate instances. Run explicit host concurrency 1,2,4 while internal filter/FFT threads stay 1. Record actual in-flight request limits. If a cache exists, compare disabled/cold/warm/evicted operation and stress retained-entry lifetime at its configured memory bound. Every fixed-configuration repeated result must be bit-identical.

Build configurations run sequentially with total compiler/test concurrency <=4, including MSVC /MP. Do not force unsupported ISA paths. Report sanitizer/platform/runtime coverage separately.

Store manifests, frozen budgets, metrics, failure/skip inventory and a human report under docs/validation/phase-2/. Link results to case IDs and record scope. Close only after both public temporal matrices, actual Highway execution, independent operator/safety checks and phase-1 regression pass against pinned references. Source inspection and mathematical examples alone do not close phase 2.
