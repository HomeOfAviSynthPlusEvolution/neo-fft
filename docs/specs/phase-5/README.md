# Phase 5: DFTTest temporal block overlap-add

Status: implementation specification, 2026-09-22; not an implementation or runtime acceptance report. Extends [phase 4](../phase-4/README.md). This batch absorbs the original AviSynth DFTTest temporal-block algorithm into the existing Neo-FFT VS CPU path. AVS entry points and final distribution move to phase 6.

| Contract | Contents |
| --- | --- |
| [Interface](dfttest/plugin.md) | tmode=1, even/odd tbsize, tosize, precedence and validation |
| [Temporal geometry](dfttest/kernel-temporal-geometry.md) | Absolute block lattice, covering blocks, boundary slots and examples |
| [Window and reconstruction](dfttest/kernel-temporal-ola.md) | Temporal normalization, PSD/mean/profile calibration, ordered synthesis |
| [Execution](execution.md) | Request closure, optional cache, deterministic scheduling, resource bounds |
| [Acceptance](acceptance.md) | Mathematical oracles, original AVS reference comparison and regression gates |

No FFT3D algorithm changes. Keep the existing DS2 bridge, PocketFFT backend, own scalar/Highway kernels, public defaults, array API, format domain, plane selection and phase-4 dither. This is not a wholesale port of the original plugin's parameter parser, FFTW loader, global working buffers or thread pool. GPU, packed YUY2, stacked lsb input/output, file-based profiles and unrestricted transform lengths are not opened.

## Fixed original reference

Use [pinterf/dfttest at 465d5d184b1be47244b7f8b4f3da7a0d57aa4e29](https://github.com/pinterf/dfttest/tree/465d5d184b1be47244b7f8b4f3da7a0d57aa4e29), the original tritical lineage maintained for AviSynth+. This immutable source is the new tmode=1 reference. Prior modes retain their fixed neo-dfttest baseline; do not change them to match a different upstream default.

| Source anchor at the fixed revision | Extracted behavior |
| --- | --- |
| dfttest/dfttest.cpp: GetFrame_T, around lines 880-949 | Negative grid anchor, block stride, endpoint mapping, overlap accumulation |
| dfttest/dfttest.cpp: mapn, around line 1514 | Clamp temporal source indices to first/last frame |
| dfttest/dfttest.cpp: constructor, around lines 1682-1700 | Mode/parity/overlap constraints |
| dfttest/dfttest.cpp: normalizeForOverlapAdd/createWindow, around lines 2578-2614 | Squared-window normalization and stored 3D window |
| dfttest/dfttest.cpp: temporal worker, around lines 575-660 | Per-volume FFT/mean/filter and multi-slice synthesis |
| dfttest/dfttest.cpp: noise estimation, around line 2167 | Unnormalized sample window distinct from filtering window |
| dfttest/dfttest.cpp header and LICENSE | GPL source notices; preserve applicable notices when adapting source |

Source was inspected at that revision; the original binary has not been built or accepted by this specification work. Implementation must record its reference build/backend/hash. Keep source attribution with adapted code and existing licensing metadata; no additional development diary is required.

## Explicit decisions

- Absorb the original absolute time-block lattice, endpoint duplication and squared-overlap normalization. Allow even T only in tmode=1.
- Retain this project's admitted T<=15 and T<=clip length, despite the original having no 15-frame upper bound. Larger transforms are outside this batch.
- Use independent output requests with bounded working memory, rather than reproducing a mutable last-block/output-ring cache. No Kalman-like history replay or shared current output state.
- Process only blocks contributing to n. The reference may compute extra preceding blocks to warm its output ring; their results cannot affect n and are not required dependencies here.
- Preserve phase-3 curve/sample-table semantics and phase-4 deterministic dither, including deliberate reference differences. The original's default sigma/tbsize, string/file inputs, packed formats, RNG and dispatch numbers are not imported.
- Preserve canonical summation order. Caching, host concurrency and reserved threads values must not change output. Scalar versus Highway/backend comparisons use calibrated budgets; execution-order comparisons within one arithmetic configuration are bitwise.

Implement scalar geometry/window/operator first, then the formal VS path and original-reference comparison, then Highway and full regression. The original AVS plugin is a test oracle in this phase; that does not require exposing Neo-FFT's AVS entry early.
