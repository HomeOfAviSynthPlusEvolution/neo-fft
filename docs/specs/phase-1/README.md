# Phase 1: shared FFT and spatial filters

Status: specification draft for implementation. Source extraction and mathematical examples are available; no implementation or VS differential acceptance has passed yet.

## Deliverable

One VS DLL registers `neo_fft.FFT3D` and `neo_fft.DFTTest`. Both functions execute complete spatial pipelines, using the same checked storage, normalized FFT adapter and scalar/Highway dispatch infrastructure. Phase 1 closes only after both functions pass the [acceptance procedure](acceptance.md) through their public VS interfaces.

| Function | Implemented subset | Required explicit argument |
| --- | --- | --- |
| FFT3D | Spatial Wiener; constant sigma; degrid; all three spatial windows; selected-plane processing | `bt=1` |
| DFTTest | Spatial ftype 0..4; zmean; both spatial modes; all twelve windows; f0beta; output conversion without dither | `tbsize=1` |

Keep the historical temporal defaults of 3. A call omitting `bt` or `tbsize` therefore fails with a phase-1 unsupported-mode error. Changing defaults to 1 would conceal a script behavior change when temporal processing becomes available.

Phase 1 targets Windows x64 and fixed-format planar GRAY, YUV and RGB with 1 or 3 planes, integer 8/10/12/14/16-bit or float32. Use actual host plane dimensions, including subsampled chroma. Alpha, packed formats, variable format/size, AVS, GPU, temporal modes, noise estimation, curves, enhancement, ROI, interlacing and dither are outside this batch. Validate unsupported requests at creation, including calls that select no planes.

## Specification set

| ID | Contract | Consumers |
| --- | --- | --- |
| MEM-001 | [Checked span2d views and geometry](kernel-plane-geometry.md) | Every host and kernel boundary |
| FFT-001 | [FFT shape, layout and normalization](kernel-fft.md) | Both pipelines and backend tests |
| WIN-001 | [Window application and reconstruction](kernel-window-reconstruction.md) | Shared operations, not shared algorithm calibration |
| RUN-001 | [Execution and errors](execution.md) | Plans, workspace, VS adapter |
| F3D-GEO-001 | [FFT3D grid and reflection](fft3d/kernel-geometry.md) | Allocation, gathering, cropping |
| F3D-WIN-001 | [FFT3D windows, conversion and reconstruction](fft3d/kernel-window-reconstruction.md) | Spatial pipeline |
| F3D-WIENER-001 | [FFT3D Wiener and degrid](fft3d/kernel-wiener.md) | Scalar/Highway spectral kernels |
| F3D-VS-001 | [FFT3D VS interface](fft3d/plugin.md) | Registration and black-box tests |
| DFT-GEO-001 | [DFTTest geometry](dfttest/kernel-geometry.md) | Both spatial modes |
| DFT-WIN-001 | [DFTTest windows and reconstruction](dfttest/kernel-window-reconstruction.md) | Analysis, calibration, synthesis, quantization |
| DFT-FILTER-001 | [DFTTest mean removal and filters](dfttest/kernel-filter.md) | ftype 0..4 |
| DFT-VS-001 | [DFTTest VS interface](dfttest/plugin.md) | Registration and black-box tests |
| AMEND-001 | [Validation hoisting and precondition guarantees](amendment-validation-hoisting.md) | Boundary verification, RealFFT, execution pipeline |
| P1-ACCEPT-001 | [VS differential acceptance](acceptance.md) | Batch closure |

These are specification boundaries, not a requirement to create one C++ file or one development batch per row. Organize specifications by phase and function, using `plugin.md` / `kernel-*.md` to define inputs/outputs, equations, domains, errors and examples. Later phases add their own function specifications and reference existing shared operators rather than duplicate them. Do not create future phase directories before their specifications are written.

## Implementation order within the batch

1. Checked views, owned buffers, direct-DFT test oracle and PocketFFT adapter.
2. Both geometries and window/reconstruction paths; confirm their different amplitude scales.
3. Scalar spectral filters, complete frame pipelines and both real VS registrations.
4. Run old/new VS comparisons; resolve specification/reference discrepancies.
5. Add Highway kernels, repeat differential tests and complete the report.

No intermediate step closes this batch. The optional FFTW adapter can be added when needed to isolate backend errors; PocketFFT is sufficient for the required same-backend comparison because both fixed references support it.

## Common numerical and memory contracts

Use `span2d::Plane<T>` directly for non-owning 2D views, with checks in the implementation's validation layer and DS2 adaptation at the bridge. This layer is to be implemented; its required checks are defined in [MEM-001](kernel-plane-geometry.md). Do not create another PlaneView type. Owned buffers and outer block/time/spectrum descriptors have separate responsibilities.

The phase-1 format domain is fixed planar GRAY/YUV/RGB, 1 or 3 planes, UInt8/10/12/14/16 and Float32. Unsupported depths/formats are errors; there is no implicit conversion. Unselected samples are copied bitwise. Each algorithm defines its own conversion, windows, power calibration and final clipping.

Working samples, spectra and accumulators use binary32. Configuration/window derivation uses the precision explicitly specified by each algorithm: FFT3D window generation remains float; DFTTest window generation uses double before storing float. Baseline round-to-nearest-even arithmetic does not imply ties-to-even pixel quantization. Disable global fast-math and implicit contraction; later approximations need scoped evidence.

Forward FFT has scale 1 and inverse FFT scale 1/N. Batch count is excluded from N. DFTTest's legacy window scaling requires an explicit synthesis factor N; FFT3D instead removes its legacy post-inverse 1/N. Neither algorithm applies a generic run-time division by overlap weight.

Reject non-finite parameters at creation and non-finite selected samples/intermediates during frame evaluation. Preserve specified finite out-of-range behavior: FFT3D float output clips [0,1], while DFTTest does not clip. Never rely on an out-of-range float-to-integer cast for saturation.

Every published frame is fully initialized and inherits the current source frame's properties. Stateless outputs are independent of request order and concurrency within a fixed dispatch/build. Error text owns its storage across callback return; failure releases resources and never publishes a partially successful frame.

See [phase-1 acceptance](acceptance.md) for the required independent math and old/new VS comparisons.

## Reference versions

Behavior was extracted from the following committed sources. The local reference filters reside under DualSynth2's `external_projects/` directory. Their unpublished Highway modernization is not itself evidence of runtime equivalence.

| Reference repository | Commit |
| --- | --- |
| `HomeOfAviSynthPlusEvolution/neo_FFT3D` | `bc9dd39410d06470d55ffd9147a98b9252f451f3` |
| `HomeOfAviSynthPlusEvolution/neo_DFTTest` | `e5d59aaff79aa45a8dc7ccb5501aff39dd489475` |
| `HomeOfAviSynthPlusEvolution/dualsynth2` | `2a24d6b4fe808692bfa10c1f9734a3c50c15774e` |

Use these revisions for reference builds and record their binary hashes and runtime settings in the acceptance results. The equations, input domains and errors in the operator/interface specifications define the new implementation; reference behavior outside those domains is not an additional requirement. Do not copy or link reference kernels into the new implementation or use them as the sole mathematical oracle.
