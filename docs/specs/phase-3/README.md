# Phase 3: noise models and spectral enhancement

Status: implementation specification, 2026-09-22; not a runtime acceptance report. Extends [phase 2](../phase-2/README.md), superseding only the restrictions explicitly opened here. Parameter registration order, formats, geometry, temporal parity/boundaries, FFT normalization, reconstruction, properties and unselected-plane copies remain unchanged.

| Contract | Contents |
| --- | --- |
| [FFT3D interface](fft3d/plugin.md) | sigma2..4, sampled patterns, pshow, enhancement, bt=-1 |
| [FFT3D noise](fft3d/kernel-noise.md) | Frequency coordinates, model selection, automatic/manual samples, power calibration |
| [FFT3D enhancement](fft3d/kernel-enhancement.md) | Sharpen/dehalo equations and placement relative to Wiener/degrid |
| [DFTTest interface](dfttest/plugin.md) | Curves, ssystem, nlocation, alpha |
| [DFTTest profiles](dfttest/kernel-profile.md) | Interpolation, independent noise window, zmean, energy conversion |
| [Execution](execution.md) | Dependencies, immutable publication, failures, memory bounds |
| [Acceptance](acceptance.md) | Mathematical oracles, reference comparisons, regression matrix |

Implement scalar tables/operators and both public VS paths first, then Highway with scalar comparison. No AVS, Kalman (bt=0), ROI, interlacing, dither, GPU or DFTTest tmode=1 is opened. FFT3D bt=2 and 4 remain valid; DFTTest T remains odd, 1..15 and <=N.

## Reference and explicit decisions

Source extraction uses the commits pinned in [reference sources](../../reference/sources.md): FFT3D bc9dd39410d06470d55ffd9147a98b9252f451f3 and DFTTest e5d59aaff79aa45a8dc7ccb5501aff39dd489475. Both local HEADs were verified on the date above. Paths below are relative to those repositories; untracked reference material does not define behavior.

| Contract | Pinned source anchors |
| --- | --- |
| FFT3D profile | src/engine/pattern_analysis.cpp: SigmasToPattern, FindPatternBlock, SetPattern; src/engine/fft3d_engine.cpp: initialize_pattern_buffers, handle_pattern_frame |
| FFT3D enhancement | src/engine/fft3d_engine.cpp: initialize_sharpen_window, initialize_dehalo_window, normalize_params_for_sample_format, process_wiener_frame; src/code_impl/Apply.cpp, Sharpen.cpp |
| DFTTest curves | src/core_window.cpp: parseSigmaLocation, interp, getSVal; src/engine/dfttest_engine.cpp: initialize_spatially_varying_sigmas |
| DFTTest samples | src/engine/dfttest_engine.cpp: prepare_noise_points, build_noise_profile_if_needed; src/core_window.cpp: createWindow |

The following are project decisions, not claims of reference equivalence:

- Positive FFT3D pfactor multiplies sampled **power** linearly. The pinned source only tests zero/nonzero and never multiplies the pattern by its magnitude. Reference equality is required at pfactor=1; other strengths use the independent oracle specified here.
- pshow uses immutable analysis windows and private rectangular preview windows. The pinned implementation mutates engine windows; repeated output must instead be request-order-independent. Preview uses the selected plane's normal chroma flag; the reference's second transform uses non-RGB as a flag, also affecting luma.
- Automatic search with no interior candidate is rejected; valid candidates are considered even above the reference's initial 1e15 score sentinel. Out-of-grid manual coordinates are rejected.
- Curve positions must be unique after binary32 conversion. All-singleton DFTTest dimensions use an explicit DC rule instead of division by zero. Duplicate-knot and zero-dimensional reference corners are not compatibility targets.
- Validate raw domains, finite conversion, array structure and static DFTTest sample rectangles even for inactive arrays. Derived arithmetic and active FFT3D sample-grid validation follow the [execution matrix](execution.md#validation-before-derived-table-construction). Do not reproduce unsafe indexing, mutable-table races or uninitialized preview coordinates.

Preserve finite historical arithmetic quirks where specified: FFT3D has three different frequency grids, its ht is not bit-depth-scaled, and DFTTest ssystem=1 with only axis arrays uses the temporal curve. These are tested contracts, not interchangeable generic radial profiles.
