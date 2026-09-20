# Phase-1 acceptance through VapourSynth

Specification: P1-ACCEPT-001. Status: acceptance protocol, not test results. Numeric budgets, binary identities and runtime reports remain to be produced during implementation. [Phase 1](README.md) stays open until every required gate below is met.

## Comparison runs

Use clean processes with explicit DLL loading and recorded VS/Python versions. Resolve and record the actual legacy namespace/function registrations from each reference DLL; do not guess names from filenames. Never load an unrecorded autoload plugin as the oracle. Execute the **public VS function** and request resulting frames; core-only tests do not satisfy this gate.

For every mandatory compatible case, run at least:

| Run | Reference | New build | Purpose |
| --- | --- | --- | --- |
| A | Pinned neo, PocketFFT, opt=1 | PocketFFT, opt=1 | Scalar algorithm and normalization comparison |
| B | Same reference as A | PocketFFT, opt=0 | New Highway versus fixed old output |
| C | Pinned neo, PocketFFT, opt=0 | New PocketFFT, opt=0 | Compare to the actual modernized reference dispatch |
| D | New opt=1 | New opt=0 | Attribute optimization differences separately |

Set old and new internal worker/FFT threads to 1 where supported, FFT3D mt=false, dither=0 and explicit spatial temporal-size=1. Capture effective backend, CPU dispatch and planning controls. If A and C disagree, retain both outputs and resolve the reference discrepancy; do not quietly choose the easier oracle. On a host without an available non-scalar Highway target, report the missing SIMD run and leave that required validation pending.

Optional FFTW runs isolate backend error and do not replace required same-PocketFFT comparisons. Running an earlier released plugin is supplemental and must not be labeled as the pinned modernized reference.

## Reproducible inputs

Generate source samples directly in VS frames so external resize/colorspace filters cannot change them. Save the generator version and sample hashes. Include zero, constants, integer chroma midpoint, signed float chroma, corner/center impulses, horizontal and vertical ramps, checkerboards, horizontal/vertical single frequencies and deterministic noise. Include spatially varying content plus different constants/noise per frame to detect stale workspaces.

For noise use uint32 LCG state `s=(1664525*s+1013904223) mod 2^32` and derive samples from its upper bits; define row-major plane/frame traversal in the fixture manifest. Calibration uses seed 1; acceptance holdout uses seeds 17 and 101. Never tune error budgets using holdout outputs. Unit tests can use independent analytic values and do not need the same generator.

At minimum use:

- GRAY, YUV444, YUV422, YUV420 and RGB, each at 8/10/12/14/16-bit integer and float32, for the baseline spatial configuration and plane-copy checks. Construct/query exact supported VS formats; unavailable formats are explicit missing coverage, not passing skips.
- Main dimensions 128x96 for GRAY/RGB/YUV444 and 256x192 for subsampled formats. Also include a non-divisible size such as 130x98 (260x196 for subsampled input), and kernel-specific first admitted and first rejected small sizes.
- Clips of 1 and 7 frames. Request all frames sequentially, then shuffled and with repeated indices. Run separate instances and concurrent requests in the new plugin; repeated outputs within one fixed dispatch configuration must match exactly.
- Float inputs including -0.25,0,0.5,1,1.25 and selected-plane NaN/Inf rejection. Use copied float planes with unusual bit patterns to verify copying does not filter or canonicalize them.

## Functional matrix

Each row is mandatory. Use pairwise parameter combinations plus explicitly listed interactions, instead of claiming an exhaustive Cartesian product. Every selectable feature must occur in a real VS comparison. Record actual case IDs and resolved parameters, not only a total test count.

| Case family | Required configurations | Observable checks |
| --- | --- | --- |
| P1-F3D-BASE | bt=1; sigma 0,2,8; beta 1,2; degrid 0,0.5,1 | Both algorithm defaults and nontrivial Wiener attenuation; constant restoration; no ignored degrid |
| P1-F3D-GRID | Blocks 8x8,8x6,9x7,32x32; overlaps zero, negative default, legal half and differing X/Y | Border/crop alignment, odd/even FFT, every visible edge/corner |
| P1-F3D-WIN | wintype 0,1,2 with nonzero overlap and degrid on/off | Window/filter/reconstruction interactions, not just identity |
| P1-DFT-FILTER | ftype 0..4, zmean false/true, both smode values | Every filter through VS; sigma2 effect for type 3; pmin/pmax bands for 3/4 |
| P1-DFT-VALUES | ftype 0/1 sigma 0,8,32; type 2 sigma 0,0.5,1; type 3 sigma 0.25/sigma2 1.5 and distinct bands; type 4 nonzero sigma,pmax | Distinguish thresholds from gain units; avoid testing only identity or all-zero output |
| P1-DFT-BETA | ftype=0; f0beta 0.5,1,2 and representable neighbors of both branch boundaries | General power, sqrt and direct branch selection; core tests also construct exact powers |
| P1-DFT-GRID | smode=1 B=8,9,16; O=0,floor(B/2), and B=8/O=6; smode=0 B=1,3,9 | Effective overlap normalization, heavy overlap, center-sample vs overlap-add |
| P1-DFT-WIN | swin 0..11; twin 0..11 at T=1; Kaiser beta 0 and 2.5; both smode values | Full window formulas, especially asymmetric ID 9 and non-unit center gain |
| P1-PLANES | Omitted, empty, [0], [1] where present, all indices and duplicates | Different empty-list rules; bitwise unselected samples and unchanged properties |
| P1-ERRORS | Invalid type/index/enum, overflow/non-finite parameters, invalid overlap, unsupported modes and execution controls | Creation vs frame error; no fake success/pass-through; restrictions classified separately |
| P1-META | Sentinel user properties plus color/range/timing properties | Same visible dimensions/format/count/rate and preserved properties |

Cover nontrivial FFT3D window/degrid and DFTTest window/zmean/ftype interactions on Gray float and high-bit YUV in addition to the all-format baseline. For large/high-overlap cases the fixture may use a small admitted image to bound test time, but it must still hit all requested branches.

Respect each equation's admitted domain when generating combinations. For example flat-top B=4/O=0 with zmean=true has zero template DC and belongs to the creation-error matrix, while the same window with zmean=false belongs to numerical comparison. Invalid combinations cannot be counted as successful filtering coverage.

## Independent lower-level checks

VS comparison is necessary but not sufficient. Also verify FFT sign/shape/scaling against binary64 direct DFT; the documented hand examples; coverage weights; units/mean restoration; exact spectral threshold equality and adjacent values; scalar/Highway tails with guard or sanitizer evidence; checked geometry overflow; natural-but-not-vector alignment; row padding; active batch counts 0,1,C-1,C; failure cleanup and reuse.

Do not force unsupported CPU instructions. Record actual Highway target, supported ISA list and sanitizers/platforms exercised. Build configurations run sequentially with total compiler concurrency <=4. A successful cross-compile does not establish runtime safety.

## Numerical acceptance

Exact gates: copied visible samples, structural metadata and property payloads, error stage/category on the agreed domain, and new-build repeatability. Padding bytes are not part of frame equality and must not be read to calculate a hash.

For processed integer planes, record maximum absolute LSB difference, differing-pixel fraction, MAE, RMSE and locations of extrema. For float, record those applicable metrics plus per-sample `abs(a-b) <= atol + rtol*max(abs(a),abs(b))`. Report each plane and separate full-frame/interior/edge-band/corner results. Edge-band width is the larger relevant block dimension, capped to the image; full-plane maxima are always retained. PSNR alone is not a gate.

Before examining acceptance holdout results, produce and freeze a budget file keyed by algorithm, sample format/units, backend pair and dispatch pair. Each processed-output comparison needs explicit finite max-error limits (integer LSB, or float atol/rtol) and any fraction/RMSE limits used. Calibrate against direct-DFT error and reference/new scalar diagnostics on the calibration fixtures, explaining the numerical source of each allowance. Avoid arbitrary universal 1-LSB/PSNR promises and avoid unlimited relative error near zero. **An absent or unreviewed budget blocks closure.**

Threshold discontinuities require branch-sensitive core cases in addition to pixel budgets. For end-to-end threshold-sensitive failures, inspect the corresponding residual/power and reference path: systematic wrong scale, wrong branch, wrong border or wrong window cannot be waived by loosening a global tolerance. Changing a frozen budget requires documented numerical evidence, a versioned decision and a fresh independent holdout, not simply rerunning until green.

## Report and closure

Store a machine-readable manifest and concise human report under a future `docs/validation/phase-1/` directory. Each result links its spec ID and includes:

1. Old/new source commits and dirty status, build recipe/toolchain/dependency revisions, DLL paths and SHA-256, VS/Python versions, CPU/OS/ISA and effective thread/backend settings.
2. Input generator/version/seed/hash, function/registration, explicit parameters, output frame/plane hashes, metadata/error captures and raw metric results.
3. Frozen budget version, calibration/holdout separation, failed or skipped cases, explained compatibility exceptions and remaining scope.
4. Separate conclusions for mathematical checks, scalar compatibility, Highway compatibility, concurrency/memory and optional platforms/backends.

Close only when both public functions' required spatial matrices pass old/new VS comparisons, Highway checks pass on an actually exercised target, required core/safety checks pass, and every discrepancy is resolved against the corresponding operator/interface specification and documented in the test results. A missing reference DLL/runtime, skipped VS tests or only new-scalar/new-Highway agreement leaves phase 1 open. AVS absence does not block phase 1; AVS matching belongs to the final batch.
