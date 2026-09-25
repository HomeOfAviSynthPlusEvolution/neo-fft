# Execution, resources and precision boundaries

[Contents](../README.md)

Each output request consumes declared input frames and owns mutable workspace. FFT3D planes and DFTTest blocks run on the host calling thread; PocketFFT is also internally single-threaded. There is no internal thread pool. The host may issue concurrent output requests; legacy execution arguments such as threads exist only in host signatures, are never read or forwarded to the core, and affect neither scheduling nor retention.

## What makes results deterministic

Configuration, windows, FFT plans and published noise models are immutable. Cached input spectra are read-only; Kalman checkpoints are restored into request-private state. DFTTest accumulates in temporal-start, spatial-Y, spatial-X order; FFT3D preserves horizontal-then-vertical synthesis. Dither uses output coordinates rather than worker identity.

Within one build/backend/dispatch/arithmetic configuration, non-Kalman modes should agree under sequential/concurrent requests and cold/warm caches. Kalman is an exception: kalman_warmup bounds cold work, while available checkpoints and request history select the starting state. Eviction or concurrent scheduling can change output. Private mutable state provides ownership safety, not request-order-independent pixels.

Changing scalar/Highway kernels, FFT profiles or arithmetic grouping can change rounding. Evaluate integer LSB, float absolute error, threshold decisions and clipping separately; there is no universal project-wide epsilon. Reference comparisons also depend on modes, sample units and random-noise definitions. Historical test budgets are evidence for those cases, not a theorem for all inputs.

## Three distinct budgets

| Object | Retention rule | Excludes |
|---|---|---|
| Idle workspace | At most one across all planes per instance, at most 64 MiB | Independent active workspaces and host frames |
| FFT3D raw spectrum rows | Default cache_mb 128 MiB, also limited by cache_frames | Workspace, host cache, request registrations and other exclusions in the cache article |
| Kalman checkpoints | Default max(64 MiB, one complete checkpoint plus metadata), at most four | Evicted snapshots still held by requests, and request-private state |

Four concurrent host requests can allocate four active workspaces despite an idle-retention count of one. That limit does not serialize them or cap process memory at 64 MiB. Spectrum-cache admission failure uses private computation; actual allocation failure is still an error.

## Lifetime and failure

Read source frames only while holding valid host ownership; do not fetch under cache locks. Models and spectrum rows publish only when complete. Failed outputs and partial Kalman states are not published. Staged Kalman replay releases the previous source frame, but host history metadata and upstream caches are outside plugin memory guarantees.

Type/domain/geometry/derived-table errors occur at creation. Active sample, FFT or output arithmetic overflow fails frame evaluation. Unselected planes are copies except explicit AVS mode 1, which skips writes; algorithm parameters still undergo applicable raw validation when inactive, except completely ignored execution compatibility arguments. Output properties always come from source n.

See [API execution rules](../../../api/en/README.md), [spectrum cache](../fft3d/spectra-cache.md), [Kalman](../fft3d/kalman.md) and [KernelInfo](../kernel-info.md).
