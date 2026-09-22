# Phase-3 dependencies, model lifetime and kernels

Specification: RUN-003. Extends [phase-2 execution](../phase-2/execution.md). Estimation is source-derived initialization, not output-history recurrence.

## Request closure

Compute an immutable request description with effective modes and logical slots. Union real indices for host requests, retaining duplicate logical transform/sample slots.

| Mode | Declared indices |
| --- | --- |
| FFT3D uniform/analytic denoising | Phase-2 neighborhood/fallback for n |
| FFT3D sampled denoising | That neighborhood union {clamped pframe} |
| FFT3D bt=-1 without preview | {n} |
| FFT3D effective pshow | {n}, irrespective of bt/pframe |
| DFTTest ordinary/curved/inactive nlocation | Phase-2 clamped T slots |
| DFTTest active nlocation | Those slots union every [fn,fn+T-1] sample interval |
| DFTTest planes=[] | {n} only |

Retain this declared closure even when a model cache is warm, so request/process agreement cannot race publication. Warm models need not reread sample pixels. No undeclared fetch is allowed. Read scope comprises output-selected planes plus named DFTTest sample patches, or FFT3D selected sampling planes. Non-finite values outside consumed regions/planes must not fail a frame merely because the containing frame was requested. Host failure delivering a declared dependency may still fail the request.

No synchronous frame fetch at creation. Validate static indices, rectangles and configuration arithmetic there; gather models after host dependencies arrive. Properties always come from n, even if deduplicated with a sample frame. Use checked frame arithmetic without adding a new arbitrary clip-length cap in this phase.

## Validation before derived-table construction

First validate all supplied types, finite conversion, scalar domains and array structure, including overridden controls. DFTTest sample tuple indices/rectangles and curve endpoints/unique positions are always checked, including planes=[] and ftype>=2. FFT3D manual bounds or an automatic candidate set are required only when sampling/preview actually runs. An analytic model or bt=-1 without preview does not activate that search.

Resolve the effective mode next, then construct only its consumed tables/constants. The finite scaled-sigma values needed for FFT3D's zero-pfactor model selection are still checked; positive pfactor selects sampled mode without computing unused scaled sigmas. Inactive *derived* arithmetic does not create an error: this explicitly overrides any inherited requirement to eagerly create all possible tables.

| Effective path | Derived work required |
| --- | --- |
| FFT3D effective preview | Normal analysis window/grid and pcutoff weights for selection, private rectangular preview windows; no noise-strength table or enhancement constants |
| FFT3D bt=-1 without preview | Normal reconstruction; only enabled enhancement constants/windows and its grid when degrid is enabled; no Wiener noise table, sample weights or sample fetch |
| FFT3D denoising | Only selected uniform/analytic/sampled model; only enabled enhancement gains; sampled model ignores sigma1..4, analytic/uniform models ignore sample weights |
| DFTTest active nlocation | Normal filtering window/template plus h2/G2 and estimated A; validate raw curve arrays but do not raise/interpolate them or derive an overridden scalar A |
| DFTTest without active nlocation | Normal filtering and active scalar/curve A; no h2/G2 or estimated table |
| DFTTest planes=[] | Raw configuration/sample/curve validation and checked copies; no FFT/window/profile allocation |

Within an active gain, reject non-finite derived constants even if a particular frame might avoid their use. Source-dependent failures, including non-finite sampled power or estimated A, occur during frame evaluation. Static window/template failures occur at creation. Thus a finite large ht whose square overflows is admitted when dehalo=0 and rejected at creation when dehalo>0; a structurally valid curve whose derived values overflow is irrelevant when active nlocation replaces it. Finite raw values outside their declared domains remain invalid even when unused.

## Publication and bounds

Analytic/curve tables are immutable at creation. Sampled models have a synchronized absent/building/ready lifecycle. Build privately, validate completely, then publish a shared immutable object with a happens-before relation. Never use an unsynchronized ready flag or overwrite live sigma arrays. Deterministic duplicate builds are permissible; only a complete model may win publication.

On failure discard partial models, release leases and fail that request; leave initialization retryable. Never reuse partial sums or poison unrelated requests. Do not hold model locks while fetching/waiting for host frames. Destruction respects active model/request leases.

FFT3D stores one base P per selected plane, with request-specific T*P. DFTTest stores one shared T*S*K table across selected output planes. Models belong to a source/filter/configuration, never only a frame index. Preview needs no persistent model. Do not allocate clip-length-proportional caches or retain all sample host frames. Gather tuples sequentially with scratch bounded by a transform, window/template, accumulator table and existing workspace, not M*T full-frame copies. Declared lists can contain up to 500*T sample indices before deduplication.

## Kernel boundaries

Separate model generation from filtering. Pass checked table views or explicit uniform descriptors, bin count/shape and aliasing contracts. Temporal FFT3D noise is spatial and broadcast to time; DFTTest tables cover the full 3D half-spectrum. Do not infer the model from ambiguous null pointers or process padding as bins. Capture block DC before writes. Odd-dimension SIMD tails cannot depend on hidden padding.

Spectra, enhancement intermediates, removed means, preview windows and inverse/accumulation buffers are request-private. Pool budgets include new simultaneously live buffers with checked byte arithmetic. Concurrent/random-order calls remain supported. Own kernels are scalar for opt=1; this does not force a scalar FFT backend. No kernel calls host APIs or starts workers.
