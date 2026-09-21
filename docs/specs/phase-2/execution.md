# Scheduling, dependencies and temporal workspace

Specification: RUN-TIME-001. Extends [RUN-001](../phase-1/execution.md) and the checked memory/ownership contracts. Let N be positive clip frame count and n the requested index in [0,N-1]. Calculate index arithmetic in a checked wider type before narrowing to the host API.

## FFT3D slots and exact fallback

For configured bt in 1..5, use integer divisions:

```text
left = bt/2
right = (bt-1)/2
T = 1 if (n<left or N-1-n<right) else bt
c = T/2
slot[j] = n-c+j, j=0..T-1
```

The boundary decision precedes frame requests, noise selection and cache lookup. Do not request out-of-range frames, clamp neighbors, choose another partial temporal length or mutate the instance's configured bt. T=1 needs only clip[n].

| Configured bt | Full window | Frames using temporal filtering |
| --- | --- | --- |
| 2 | n-1,n | 1<=n<=N-1 |
| 3 | n-1,n,n+1 | 1<=n<=N-2 |
| 4 | n-2,n-1,n,n+1 | 2<=n<=N-2 |
| 5 | n-2,n-1,n,n+1,n+2 | 2<=n<=N-3 |

An empty interval means every frame falls back. With N<bt, all frames use 2D. In particular the last frame for bt=2 can use temporal filtering: 'first and last frames always fall back' is incorrect.

## DFTTest slots and endpoint duplication

Validate odd T=tbsize in 1..15 and T<=N at creation, independently of n. Set c=T/2 and retain exactly T slots:

```text
nominal[j] = n-c+j
slot[j] = clamp(nominal[j],0,N-1), j=0..T-1
requested_frames = unique(slot[0..T-1])
```

The host dependency set may be deduplicated, but numerical slots **must not** be deduplicated or compacted. Their length and ordering stay T. Center slot c always refers to n, even when earlier/later slots refer to that same endpoint. Apply tw[j] by slot index, not by unique-frame index.

Examples for T=5,N>=5: n=0 gives [0,0,0,1,2]; n=1 gives [0,0,1,2,3]; n=N-1 gives [N-3,N-2,N-1,N-1,N-1]. For N=T=3,n=0 the slots are [0,0,1]. Reject N<T at creation; replication does not legalize a short clip. DFTTest planes=[] needs only n after configuration validation.

## Host lifecycle

At the host request stage declare the complete unique source dependency set. At the ready stage acquire all needed frame handles and retain them through every borrowed view's last use. Use the same slot map at both stages; never synchronously ask for an undeclared neighbor while holding a workspace/cache lock. DS2 may express these stages, but kernels receive only validated views and immutable configuration.

Initialize output from metadata/properties of clip[n], copy unselected planes from n, and publish after all selected planes succeed. Frame requests in increasing, decreasing, repeated or simultaneous order have identical results for a fixed dispatch configuration. Separate instances share no mutable numerical state. Missing upstream frames, allocation failures and FFT errors release every frame/workspace handle and return a stable owned error; they do not poison later requests.

There are no phase-2 nlocation/pattern dependencies. No filtered output frames are dependencies, and no lock enforces output-frame order.

## Workspace accounting

All sizes below are minima for dense logical storage, not permission to ignore actual row/plane/batch strides. Check products/offsets and account for backend scratch, frame references and output memory. Workspace pooling may reuse storage only under exclusive leases, reset accumulation for every request/plane, and must release leases on all error paths.

| Buffer | Minimum logical payload |
| --- | --- |
| DFTTest padded temporal source, if materialized | T * pad_height * pad_width * sample_bytes per active plane; allocate by actual row/plane pitch |
| DFTTest spatial accumulator | pad_height * pad_width * sizeof(float), **2D only** |
| DFTTest real transform batch | C * (T*B*B) floats, plus any explicit adapter padding |
| DFTTest complex transform batch | C * (T*B*(floor(B/2)+1)) complex<float> |
| DFTTest mean restoration | One g per active block plus immutable G, or a full private M per block; preserve g before overwriting DC |
| FFT3D temporal input batch | T * C * bh * (floor(bw/2)+1) complex<float>, or equivalent retained immutable cache slices |
| FFT3D output batch | C * bh * (floor(bw/2)+1) complex<float> |
| FFT3D temporal scratch | T complex values per concurrently processed spatial bin, plus any retained grid component |

DFTTest immutable h and G are one full volume/half-spectrum each. Spatial streaming is allowed if it preserves source lifetime, checked addressing and synthesis order. Avoid allocation proportional to total clip length. Both algorithms retain phase-1 cover/accumulator and backend scratch requirements in addition to the table.

## Optional bounded spatial FFT cache

A cache is an optimization, not a prerequisite or a source of numerical state. FFT3D cache entries contain raw windowed spatial spectra **before Wiener/degrid**, identified by source frame, plane, block region/layout and immutable format/window/FFT-plan configuration (an instance-local key may imply the latter). Do not mix different instances/configurations.

Publish entries only after successful completion; published payloads are immutable. Concurrent requests retain ownership until use finishes, so eviction cannot free borrowed storage. Give the cache an explicit entry/byte limit, checked accounting and an eviction policy independent of numerical results. Cache-off, warm-cache, cold-cache and eviction runs must agree. No numerical work may mutate an entry in place; a current-frame-derived grid correction is request-specific. Synchronization guards publication/lifetime only and cannot wait for host frames while holding a lock needed by another request.

Record measured workspace and cache limits in the implementation report. Shared backend plans must either be documented as concurrently executable on disjoint arrays or be leased/isolated safely; concurrent planning/destruction follows phase-1 backend rules.
