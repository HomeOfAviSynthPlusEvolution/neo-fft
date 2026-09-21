# Phase 2: stateless temporal filtering

Status: implementation specification. Algorithm definitions use the pinned versions in [reference sources](../../reference/sources.md). Runtime acceptance follows the protocol below.

Phase 2 extends the [Phase 1 contracts](../phase-1/README.md); it retains their sample units, spatial geometry, window equations, formats, planes, quantization, checked views, ownership and failures. The interface addenda override the old temporal-size restrictions. All other phase-1 parameter restrictions remain in force.

| Order | Contract | Purpose |
| --- | --- | --- |
| 1 | [FFT3D interface](fft3d/plugin.md), [DFTTest interface](dfttest/plugin.md) | Defaults, admitted domain and creation errors |
| 2 | [Scheduling and memory](execution.md) | Exact frame-slot mapping, boundaries and lifetime |
| 3 | [3D real FFT](kernel-fft.md) | Shape, strides, half-spectrum, sign and normalization |
| 4 | [FFT3D temporal operator](fft3d/kernel-temporal.md) | Current-frame degrid, temporal DFT, PSD Wiener and center inverse |
| 5 | [DFTTest temporal operator](dfttest/kernel-temporal.md) | 3D window energy, zmean, all five filters and synthesis |
| 6 | [Acceptance](acceptance.md) | Independent numerical checks, public VS comparison and stress |

Implement scalar operators and compose both public plugins before optimizing the new paths with Highway. bt=1 and tbsize=1 must regress against phase 1. Defaults bt=3 and tbsize=3 now work when the input satisfies the respective interface rules.

## Scope and invariants

- FFT3D supports bt=1..5; an unavailable full neighborhood falls back to btcur=1. Even lengths 2 and 4 are valid and past-biased.
- DFTTest supports **odd tbsize=1..15**, tbsize<=clip.num_frames, tmode=0, effective tosize=0. Endpoint replication fills exactly T slots; it never shortens the transform.
- Both return only frame n. FFT3D selectively inverts time at its current-frame index; DFTTest reconstructs the center slice of a normalized 3D inverse. Neither performs temporal overlap-add.
- DFTTest evaluates the exact phase-1 window formulas at length T. Center gain is the actual tw[T/2]^2, not universally 1 (flat-top window 6 is a counterexample), and must not be renormalized away.
- PSD means squared complex magnitude in the specified transform units. DFTTest evaluates power after 3D mean removal; FFT3D after temporal-DC grid adjustment. Their sigma calibration and epsilon placement differ.
- DFTTest nlocation/alpha noise estimation and sigma curves, and FFT3D patterns/frequency-dependent noise/enhancement, remain phase 3 per the [roadmap](../../design/roadmap.md). Ordinary per-bin PSD filtering is mandatory here.

Frame n is a pure function of immutable configuration and declared source frames. No output history, recursive accumulator or request-order-dependent numerical state is allowed. Immutable tables and optional bounded caches are permitted; synchronization for publication does not imply temporal recurrence. Mutable spectra, removed components, inverse buffers and accumulators belong to the active request. Kernels do not call host APIs or spawn workers.

Binary32 storage/arithmetic and binary64 window evaluation follow phase 1. Check dimensions, volume/stride/byte products and frame-index arithmetic before allocation/access. Non-finite configuration/table calculations fail creation; non-finite selected input/intermediates fail the frame without publishing a partial output.
