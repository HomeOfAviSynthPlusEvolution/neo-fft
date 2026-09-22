# Phase-3 acceptance

Specification: P3-ACCEPT-003. Extends [phase-2 acceptance](../phase-2/acceptance.md). These specifications do not claim goldens or runtime passes.

## Independent mathematical checks

Use direct DFT for small volumes and separately written table/formula oracles. Expected values must not call production code. Compare tables/intermediates as well as quantized images.

| Area | Required cases |
| --- | --- |
| FFT3D analytic | Equal sigmas, each unequal anchor, knots and neighbors, rectangular/odd blocks, K denominator, odd-H asymmetry, squared amplitude scaling |
| Pattern search | Known least-power interior block, first-winner ties, excluded border, all scores >1e15, smallest 5x5 grid, no candidates rejected, manual edge block, degrid 0/1, padding excluded |
| Sample power | Direct FFT, pcutoff/DC weights, pfactor 0.5/1/2 linear power, sigma inactive with positive pfactor, T rather than T^2 calibration |
| Temporal profile | bt=1..5, even neighborhoods, fallback then interior and reverse, no first-request calibration leak |
| Enhancement | Sharpen/dehalo/both/zero, DC, zero smin/smax/ht/svr, degenerate active window rejection, odd dy, scaled smin/smax versus unscaled ht |
| Composition | Uniform spatial fused versus pattern two-stage formula, temporal inverse before enhancement, second-stage DC ratio, bt=-1 independent of sigma, preview priority |
| Curves | Unsorted valid pairs, exact knots, root-before-interpolation, shared/separable/radial modes, radial axis-only quirk, singleton/all-singleton axes, odd T coordinates |
| Curve filters | ftype=0..4, type-dependent wscale, unchanged scalar sigma2/pmin/pmax, no sigma square, slocation precedence |
| Noise estimate | One/multiple/repeated tuples, start not center frame, formats and chroma/RGB, raw coordinates, full T transform, unnormalized h2, energy ratio, alpha linear power, zmean on/off and DC |
| Shared model | Sample unselected output plane, mixed-plane tuples form one average, all output planes use same table, ftype>=2 inactive samples, planes=[] never samples |
| Inactive derived values | FFT3D bt=-1, sharpen=dehalo=0, sigma=1e30 reconstructs without computing sigma squared, while otherwise equivalent bt=1 rejects overflow; valid large ht with dehalo=0 versus active overflow rejection; nlocation overrides valid curves without constructing them; bt=-1 and analytic noise do not require an automatic sampling grid; planes=[] performs no model allocation |
| Reconstruction | bt=-1 zero strengths; DFTTest unit-gain curves with inherited center-window gain; empty/default controls regress phase 2 |

Hand checks: FFT3D W=H=8,x=0,y=1 gives f=sqrt(0.5)/4, hence sigma3^2/norm. DFTTest T=1,S>1, slocation=[0,0,1,16], ssystem=0, fx=fy=0.5 gives V=4 (interpolated square-root endpoints multiplied), not 8. M=2 with residual powers 4,12, wscale2/wscale=0.5, alpha=5 gives A=20.

## Public VS comparison

Build fixed reference commits from the phase-3 README; record DLL/dependency hashes, compiler, host, backend and actual CPU target. Source inspection is not a black-box pass. Compare both formal VS functions for new scalar and Highway. Existing frozen tolerances apply only to their domains; register intermediate/table tolerances before judging results. Do not relax thresholds to hide calibration errors.

Cross models with spatial/temporal paths, boundary/interior, enhancement/degrid/zmean on/off, integer 8/10/16 and Float32, GRAY/YUV420/YUV444/RGB and plane selection. Include rectangular/odd/non-power-of-two blocks, T=3/5/15, bt=2/4, short clips and padded strides. Constant, impulse, exact-bin sinusoid, ramp, noise and mixed-detail synthetic patches are mandatory; external video is optional.

Classify intentional differences using the phase-3 README. pfactor!=1 uses the linear-strength oracle. Unsafe/no-candidate/duplicate/all-singleton cases use validation/math oracles. Compare pshow to fresh reference instances and separately verify the defined chroma correction. Verify the actually executed reference scalar/Highway branch, including analytic profiles with promoted effective pfactor, before treating output as an oracle. Do not accept unexplained mismatches. Exceptions do not waive ordinary pfactor=1 sampled filtering, analytic profiles, enhancement or DFTTest comparisons.

## Scheduling, errors and lifecycle

Instrument input callbacks for exact dependency unions, no undeclared reads, inclusive sample intervals, repeated logical slots/tuples despite deduplication, pframe clamp and no extra dependencies for preview/bt=-1/inactive samples/empty planes. Test simultaneous first calls, descending/random/repeated requests, deterministic models and output under the existing same-backend policy.

Reject malformed arrays, duplicate converted knots, missing endpoints, negative/NaN/Inf values, overflow, invalid sample planes/rectangles/intervals, >500 tuples, invalid grids and active degenerate windows. Sample-patch NaN on an output-unselected plane must fail; NaN in an unconsumed region must not be inspected. Simulate upstream sample errors and failed private builds; retries must not expose partial tables. Verify stable errors, lease release, properties and bitwise unselected copies.

Run scalar/Highway logical-tail checks, concurrency tests and all earlier admission/math/dispatch regressions. Build concurrency <=4, configurations sequential. Verify header dependency tracking before accepting incremental results: Plan layout changes must rebuild plugin allocation sites. Stale ABI objects are not algorithm evidence. Separate compile-only from executed platform coverage.

## Completion

Implementation closes only when public modes, independent oracles, reference comparisons with explicit exceptions, lifecycle/error checks and earlier regressions pass. Record matrix, tolerances, actual CPU/backend paths and skipped coverage. Markdown/source checks close specification work only, not implementation acceptance.
