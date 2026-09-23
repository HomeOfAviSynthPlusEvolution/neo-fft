# Phase-5 dependencies, cache and lifetime

Specification: RUN-005. Extends [phase-4 execution](../phase-4/execution.md). Only DFTTest tmode=1 changes its dependency geometry; FFT3D and DFTTest tmode=0 remain unchanged.

## Immutable request closure

At request n compute its absolute covering starts and each start's T ordered clamped slots. Host requests are the deduplicated union of these real indices plus n for properties/copies. Active nlocation adds every inherited sample interval fn..fn+T-1. Explicit planes=[] requests only n after raw validation.

Freeze this closure and the logical slot maps before host acquisition. Keep the same declared closure on warm and cold caches/models so request/process agreement does not depend on races. A warm model or block need not reread its pixels; no undeclared synchronous fetch is allowed. Duplicate logical slots still multiply by their own temporal-window entries.

The output-neighborhood set contains at most 2*T-1 real frames. The inherited maximum 500 sample tuples adds at most 500*T indices before deduplication; allocation is bounded by admitted T/sample count, never N or n. Sample tuples retain supplied order for accumulation. Finite-content checks cover selected output regions and explicitly named noise patches only; unused planes/regions and omitted noncontributing ring-warmup frames are not read. Host failure delivering a declared dependency can still fail a request.

Use the existing stateless DS2 request path or equivalent stages satisfying this immutable closure. No new DS2 API is needed beyond the fixed phase-4 dependency. There is no dependency on previously requested outputs, no replay from frame 0 and no need for a checkpoint. Request metadata always comes from n.

## Private work and optional cache

The required baseline can recompute each contributing volume per output request; a persistent block cache is optional, not an acceptance prerequisite. Keep windows, G, models and plans immutable; spectra, means, inverse buffers, output accumulators and dither rows are private/leased. Process bounded transform batches. Do not allocate a full output-history ring or N full-frame products.

If a block cache is implemented, key it by instance/configuration identity, selected plane, **unclamped temporal start**, spatial origin, and backend/dispatch/arithmetic shape. Source/model identity belongs to that instance; no reuse across a changed model/configuration. Cache immutable raw spectra or inverse block products sufficient to recover the canonical individual contributions, not already quantized/dithered frames or differently grouped spatial sums. An entry must describe exactly which processing stage it contains so mean removal, filtering or scaling cannot be applied twice.

Default optional cache bound: 64 MiB and at most 64 entries per instance, including actual allocated capacities and metadata; an entry larger than the budget is not cached. Use synchronized immutable publication and bounded LRU retention. Leased evicted entries remain alive, with their bytes counted separately in active-request diagnostics. Each request keeps only its current bounded batch of leases; do not lease the entire history. Cache absence, misses, eviction, duplicate builds and disabling the cache must leave output bit-identical in one arithmetic configuration. Existing idle-workspace budgets remain phase 4's independent limits.

No lock is held during host fetch, FFT, worker join or callback. Do not wait for another unfinished output request to fill a cache/model; deterministic duplicate computation is acceptable. Publish only complete finite products. On error discard private/partial products, release leases and fail the frame; successful unrelated immutable entries remain valid. Cancellation and teardown follow phase 4.

## Work order and bounds

Under the revised phase-4 execution policy, all transforms run on the calling host worker; threads is reserved. Commit results in temporal-start then spatial-Y then spatial-X order; within each block retain the existing pixel contribution order. Do not use unordered atomic additions or use a cache hit to change addition grouping. Dither runs once per final plane with the phase-4 deterministic seed identity.

Checked arithmetic covers signed starts/slot endpoints, floor division, T*B*B, frequency extents, strides, job counts and byte capacities. Only convert a clamped, range-checked real frame index to the host index type. A near-INT_MAX n must not trigger a loop from A to n. Memory/time for one request scales with admitted geometry, T, sample tuples and active concurrency, not distance from the first frame.

PocketFFT remains required and FFTW optional; original-reference execution may use FFTW without requiring a new production backend. Preserve all phase-4 thread/opt mappings. Performance testing may separately measure reuse benefits, but cannot justify changing block phase, reduction order, noise calibration or deterministic behavior.
