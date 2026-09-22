# Phase 4: state and complete VS CPU behavior

Status: implementation specification, 2026-09-22; not an implementation or acceptance report. Extends [phase 3](../phase-3/README.md), opening only the restrictions listed here. Public argument order, plugin identity, supported formats, metadata, phase-2 temporal rules and phase-3 model precedence remain in force unless explicitly overridden.

| Contract | Contents |
| --- | --- |
| [FFT3D interface](fft3d/plugin.md) | bt=0, kratio, ROI, interlaced, mode precedence |
| [Kalman operator](fft3d/kernel-kalman.md) | Initial state, motion reset, recurrence, zero-noise rule |
| [ROI and fields](fft3d/kernel-roi-fields.md) | Plane coordinates, row packing, geometry, copy boundaries |
| [DFTTest interface](dfttest/plugin.md) | Dither controls and format activation |
| [Dither operator](dfttest/kernel-dither.md) | Diffusion, deterministic random samples, conversion |
| [Replay](replay.md) | Canonical state, bounded cache, staged dependencies, DS2 release contract |
| [Execution](execution.md) | Worker controls, dispatch, allocations, error ownership and teardown |
| [Acceptance](acceptance.md) | Independent oracles, reference comparisons and integration gates |

Deliver both VS CPU functions, scalar kernels first and Highway parity where applicable. DFTTest tmode=1 is added in [phase 5](../phase-5/README.md); AVS integration and final distribution move to phase 6. GPU, tmode=1, a new noise estimator, new public checkpoint parameters and mandatory FFTW support remain outside this phase-4 batch. PocketFFT remains the required/default backend. FFT3D bt=-1..5 is supported here; DFTTest tbsize remains odd, 1..15 and <= clip length until the phase-5 mode extension.

## Fixed reference and intentional differences

Algorithm source anchors are relative to the repositories in [reference sources](../../reference/sources.md): FFT3D `bc9dd39410d06470d55ffd9147a98b9252f451f3` and DFTTest `e5d59aaff79aa45a8dc7ccb5501aff39dd489475`. The required Phase-4 DS2 revision is `f1d51bd0217d3878f995375e95c2827b4604facf` (`feat: release acquired staged video dependencies`), committed and pushed. DS2 contracts are defined by that revision's headers and tests, not an untracked working document.

| Area | Source anchors |
| --- | --- |
| Kalman | FFT3D src/code_impl/Kalman.cpp; src/core/cpu/kalman_hwy.cpp; src/engine/fft3d_engine.cpp: initialize_kalman_buffers, process_kalman_frame |
| ROI/fields | FFT3D src/engine/frame_cover.cpp; src/engine/fft3d_engine.cpp: configure_geometry; plugin selected-plane copy path |
| Dither | DFTTest src/core_scalar.cpp: dither<uint8_t>; src/engine/dfttest_config.cpp; workspace RNG initialization |
| Execution | Both engine configuration/worker paths; DS2 include/dualsynth/staged_video.hpp, include/dualsynth/vapoursynth/video_bridge.hpp, tests/unit/test_frame_services.cpp and tests/unit/test_vapoursynth_video_bridge.cpp |

Explicit project decisions:

- Kalman follows the reference's **sequential 0,1,...** result, including frame 0 pass-through and zero initial last spectrum. Arbitrary request order must reproduce that sequence, rather than the reference's request-history-dependent state. No hidden state reset at a seek or checkpoint eviction.
- Uniform zero noise uses the finite identity rule in the operator instead of 0/0. Sampled/analytic pattern noise keeps the reference's 1e-15 floor. In Kalman, positive pfactor selects sampled noise but does not scale its power; its magnitude remains active in inherited Wiener processing. This mode-specific rule preserves the reference sequential result.
- ROI margins are nonnegative and aligned for selected subsampled planes; mixed crop-size rounding and unsafe odd-height field packing are rejected. Interlaced processing is a full-height row permutation, **not two half-height transforms**.
- Dither mode 1 preserves the historical diffusion equations. Modes >=2 use a specified coordinate hash and default seed 0 instead of workspace/thread-local evolving mt19937 or random_device. Only UInt8 uses dither, as in the reference.
- Worker requests are resource limits/hints with a documented capability mapping. Automatic policy is deliberately fixed and does not follow the reference's changing host-concurrency heuristic.
- Bounded Kalman replay uses DS2 staged-store release and VS frame-context early release from the required revision above. Phase-4 integration must update Neo-FFT's dependency pin to that revision and verify the actual host's early-release behavior before declaring bounded replay accepted. Host bookkeeping/caches are separate from the dependency-pixel guarantee; see replay.md.

Implement in dependency order: integrate the fixed DS2 revision; ROI/packing and geometry; scalar Kalman and dither; replay/cache using the existing release API; execution controls; VS integration and acceptance. No further DS2 interface design is pending for this contract. Actual-host retention validation remains an implementation acceptance gate. Preserve phase-3 inactive-derived-arithmetic rules, with the explicit Kalman initial-covariance exception.
