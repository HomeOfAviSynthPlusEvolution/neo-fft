# FFT3D raw spectrum cache: rows, budgets and concurrent sharing

[Contents](../README.md)

Adjacent outputs reuse spatial FFTs of some source frames. This cache stores **unfiltered block spectra**, saving reflection sampling, analysis windowing and forward spatial FFT. Temporal filtering, enhancement, inverse FFT and output synthesis remain per-output work. See the [FFT3D pipeline](../fft3d.md) and [cache arguments](../../../api/en/fft3d.md#execution-and-spectrum-cache).

## 1. Why store raw spectra?

For bt=3, output 10 uses sources 9,10,11; output 11 uses 10,11,12. At a given plane and block position, source 10 and 11 spatial spectra can be reused. Cached X contains centered, analysis-windowed spatial FFT data, before degrid, temporal FFT, Wiener or enhancement. Current-block DC and neighborhoods differ between outputs, so filtered Y cannot be shared and cached X must remain immutable.

The cache belongs to one filter instance. With input, ROI, geometry, windows and FFT configuration fixed, (source frame, plane, block row) uniquely identifies an entry. Instances do not share entries. Only ordinary requests configured with bt=2..5 use it, including endpoint spatial fallback. bt=-1/0/1 and effective preview do not. Kalman checkpoints store recurrence state separately.

## 2. Row size

One entry contains spectra for **every horizontal block in one block row**. It is neither one pixel row, eight fixed blocks, nor a whole source frame. With Nx horizontal blocks and Bx×By block size, complex float32 takes eight bytes:

$$M_{row}=N_xB_y(\lfloor B_x/2\rfloor+1)\times8\quad\text{bytes},$$
$$M_{frame}=\sum_{p\in selected}N_{x,p}N_{y,p}B_y(\lfloor B_x/2\rfloor+1)\times8.$$

Block counts include padding blocks; image area divided by block area is insufficient. See [block geometry](../shared/block-geometry.md). Default blocks 32×32 and overlap 10×10 give 32×17×8=4352 bytes per block. For full-frame YUV420 with all three planes:

| Image | Y grid | Each chroma grid | One Y block row | Full source-frame payload |
|---|---:|---:|---:|---:|
| 1920×1080 | 89×51 | 46×27 | 378.25 KiB | 29.15 MiB |
| 3840×2160 | 177×100 | 89×51 | 752.25 KiB | 111.14 MiB |

Exact 1080p payload is (89×51+2×46×27)×4352=30,564,096 bytes; 4K is 116,537,856 bytes. u8, u16 and float input have identical spectrum sizes for identical geometry: cached spectra are always complex float32.

## 3. Two simultaneous budgets

Default cache_mb=128 is a demand-grown upper bound, not an allocation at creation. Small videos do not generate extra data to fill it. cache_frames limits distinct retained source-frame IDs; retaining just one row still occupies one frame slot.

Automatic VS frame capacity is bt + core thread count − 1, read at creation. AVS updates automatic capacity from thread cache hints, assuming one thread (bt source frames) until a valid hint arrives; explicit cache_frames is unaffected. bt=3 with 16 core threads gives 18. This estimates concurrency; it neither counts current in-flight requests nor preallocates 18 complete spectra. Both limits apply; legal values and disabling rules are in the [API](../../../api/en/fft3d.md#execution-and-spectrum-cache).

Charged bytes include admitted row payload, shared control objects, index nodes and distinct-frame index nodes, including rows being built or read. They exclude the fixed cache object, allocator overhead, temporary admission metadata, request registrations, active workspaces, input/output frames and host caches. Thus cache_mb=128 does not cap the process at 128 MiB.

Admission is checked before allocating a full row. If both limits permit, the row is generated and retained; if capacity cannot be freed, the caller immediately computes private blocks without waiting for space. **Budget pressure reduces reuse, not temporal length or filter behavior.** Actual allocation failure still reports an error.

The bt−1 source frames reused between sequential outputs give these payload estimates:

| Image | bt=3: two source frames | bt=5: four source frames |
|---|---:|---:|
| 1920×1080 YUV420 | 58.30 MiB | 116.59 MiB |
| 3840×2160 YUV420 | 222.28 MiB | 444.56 MiB |

These characterize working sets, not proven minimum budgets for zero redundant work. Incoming rows, metadata, concurrent requests and out-of-order access also consume capacity. When the reusable 4K set exceeds 128 MiB, eviction or bypass occurs; row granularity preserves a useful subset without retaining entire frames beyond the limit.

## 4. Joining concurrent construction

An admitted row can be building, ready or failed. The first requester produces it; concurrent same-key requests wait for that construction and share the immutable result. Waiting releases the cache mutex, allowing other keys to proceed. Producers use already-acquired source frames, with no host frame requests or recursive cache acquisition during construction; this avoids dependency cycles from joined work.

Each output holds the neighborhood rows needed for its current block row, then releases them. Entries held by producers, waiters or readers cannot be evicted. Shared ownership neither duplicates payload nor exempts it from charging.

Sixteen threads requesting **the same admitted key** can share one construction even if the budget is now full. Different keys, or an absent key that cannot be admitted, may compute privately. There is no separate unbudgeted, unlimited table sharing in-flight results.

Construction failure wakes current waiters and propagates the error without publishing partial rows. Future queries to a retained failed entry bypass it; after eviction a new admission may retry. A later filtering failure does not automatically invalidate already-complete raw spectra.

## 5. Demand and source order, without timed expiry

There is no TTL. Completing a request releases ownership and registration; unused ready rows can remain until capacity pressure or instance destruction.

Before acquiring dependency frames, ordinary host requests register their actual source interval and the next unfinished block row of each selected plane. Progress advances after each row; normal or exceptional completion unregisters it. This identifies rows still potentially needed by registered consumers.

Under byte pressure, choose among evictable rows:

1. Prefer rows with no remaining registered consumers.
2. Within that category, prefer earlier source-frame IDs, then plane and row order.
3. If every evictable row has registered demand, earlier demanded rows can still be evicted.

Frame-capacity pressure removes every retained row of a source frame. A frame with any held row is ineligible. Prefer a frame without remaining demand, then an earlier frame ID. This is not pure LRU: future registered demand is a soft priority; actual ownership is hard protection. Future random requests that have not arrived are naturally absent from the registry.

## 6. Why ideal miss rate is 1/bt

For steady sequential requests with complete neighborhoods, identical planes/rows, sufficient capacity and no intervening eviction, each output requests T source rows. T−1 are reusable and one is new. Count construction and budget bypass as misses:

$$miss=\frac{builds+bypasses}{requests},\qquad miss_{steady,ideal}=\frac1T.$$

Joining a successful in-progress construction counts as a hit. Cold starts, endpoint fallback, limited capacity, divergent concurrent row progress and random seeks change the ratio. bt=5 does not guarantee 20% for every workload.

Miss rate is not frame-time ratio. Hits save input preparation and spatial forward FFT only; temporal filtering, inverse FFT, synthesis and cache coordination remain. 1/bt therefore does not imply a bt-fold total speedup. Space figures above are payload calculations, not measured process peaks. Implementation: [row cache](../../../../src/runtime/spectra_cache.hpp) and [execution plan](../../../../src/algorithms/plan.cpp).
