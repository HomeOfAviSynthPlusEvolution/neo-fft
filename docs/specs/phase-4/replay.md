# Canonical Kalman replay and checkpoints

Specification: F3D-REPLAY-004. Applies only to effective FFT3D bt=0 without preview. Stateless modes retain phase-3 request closures. This document explicitly replaces the static phase-3 closure/publication policy for Kalman.

## State identity and request snapshot

State S_j has the meaning in [the operator](fft3d/kernel-kalman.md). Cache scope is one immutable filter instance, including source identity, selected planes, format, ROI, packing, windows, noise model, arithmetic/dispatch/backend configuration and kratio. Never share by frame index alone. Cache exact L,C,Q for **all selected planes** as one immutable checkpoint. Metadata records j and shape/configuration identity. No serialized cross-build reuse is required.

At the start of an n>0 request, under a short cache lock select and lease the largest available j<=n; if none, use virtual S0. Freeze this choice for the request. Later cache publication/eviction cannot change its dependencies or starting state. A request with j=n can render directly from cached L_n; it still needs source n for properties and copied pixels. Frame 0 requests only source 0, bypass models and state.

For j<n, clone the checkpoint into request-private mutable state, or initialize S0, then consume source frames **j+1,...,n in increasing order**. Do not consume source 0 except when separately needed as pframe or output 0. The recurrence has no finite temporal radius: evicting history means replay from an earlier checkpoint/S0, not a reset near n. One-frame clip output remains valid.

## Staged host protocol

Use DS2 staged video frame services: advance() declares dependencies, consumes completed stages and updates private state; process() produces the final output. No synchronous get-frame loop, host fetch while holding locks, recursive request for this filter's own output, or waiting for another unfinished output request.

Freeze whether a sampled model is ready together with request initialization. If absent and needed, stage the clamped pframe first, build all selected-plane models privately and publish only a complete validated immutable model. Release that source owner after sampling. Deterministic duplicate builds are allowed; never wait for another request to finish a model while holding a host worker. Uniform/analytic models need no source stage. All requests use equivalent immutable model values.

Then acquire one replay source index at a time. Validate its descriptor, consume every selected ROI, advance the whole private state, and release the acquired frame owner before declaring the next replay stage. Keep no borrowed frame/plane pointers after release. Drop temporary owning copies returned by provider get() as well: removing only the provider's copy does not release the source if callers retain owners. Source n can remain retained for final rendering; if the selected checkpoint already equals n, acquire n directly. If pframe equals a later replay index, reacquisition after release is valid. Deduplication applies to currently retained/pending owners, not to all indices ever seen.

Thus each request's source dependency sequence is the cold-model pframe if needed, followed by j+1..n, with n fetched separately when j=n. Host requests may revisit an index after release. The retained source set has at most one full-frame owner per request at these stage boundaries, independent of replay distance. Cancellation and exceptions release the current owner and private state. Index advancement uses a checked wider representation and terminates at n without overflowing n+1 at INT_MAX boundaries; there is no clip-length-sized index list or forced full-clip prefetch.

## Required DS2 version and release contract

Use DS2 `f1d51bd0217d3878f995375e95c2827b4604facf`, which is committed and pushed. Update `cmake/Dependencies.cmake` to this exact revision during Phase-4 integration; the previous Neo-FFT dependency pin predates dependency release and is insufficient. Do not fetch the removed historical revision or rely on a mutable local checkout. Without explicit release, even the new staged store retains acquired owners until the request ends: single-frame stages alone still retain O(n-j) source frames.

That DS2 revision implements `VideoStageContext::release_frame(input,n) -> Result<VideoFrameReleaseResult>` through the staged-specific `VideoStageFrameProvider`; the general provider is unchanged. Consume this API rather than adding another release abstraction. Required semantics are:

1. Release by input-node/frame identity removes that acquired owner from the request's lookup/retention set. It cannot release a pending/in-flight dependency. Invalid identity or use in the wrong stage is a defined error, not silent success or UB.
2. Release drops the staged store's owner, not owning snapshots already returned by get(). Those owners and their views remain valid until their owners expire; a bare view cannot outlive its last owner. Subsequent store lookup fails unless the dependency has been acquired again through a declared stage. Declaring it again must work; a permanent seen-index dedup set cannot suppress the fetch.
3. Host frame handles are released exactly once. Final metadata/copy paths must use a retained/reacquired n, never a dangling view. No process callback can assume all past-stage frames remain available.
4. Existing consumers that never release retain existing behavior. Both normal completion and host cancellation/error unwind owned dependencies safely. A request's acquire/release/advance state is not concurrently mutated without synchronization.
5. Retain and run the included DS2 lifetime/fake-host tests and add the Neo-FFT replay integration probe with outstanding-owner counters. Preserve the synchronous host adapter's existing behavior where that adapter is built; this does not open the Neo-FFT AVS feature work.

Release is permitted only inside advance(); invalid, unacquired or already released identities return InvalidArgument. If an output-origin dependency has been released, the existing Ready path reacquires it before output allocation/copy without invoking the filter's advance() again.

VS has a second ownership layer in frameCtx. The bridge must adopt the entire pending batch first, then call releaseFrameEarly once per unique native (node,n), preserving aliases across input indices. Missing early-release API is an explicit staged-path error. No frameCtx pointer is stored in an owner or used during owner destruction. An incomplete batch unwinds through request/host error cleanup. DS2 store release alone is insufficient to release host-held pixel references.

The locally inspected VS R73 source has a static early-release implementation risk: it traverses availableFrames using reqList.size(), while the scheduler clears reqList before completed-stage callbacks. Its history slots may also remain allocated after frame references are cleared. Actual runtime behavior must therefore be validated; API presence and ordinary staged compatibility tests do not prove bounded long-replay retention. Do not declare production bounded-replay support on an unverified host. Record a verified runtime/version condition, or resolve the host defect before enabling that guarantee; never silently truncate history as a workaround.

A local unpinned dependency override or direct raw-host bypass is not a release solution. Acceptance must record the required immutable DS2 revision and confirm the project actually builds against it. DS2 unit/fake-host and ordinary staged integration results do not substitute for the actual-host long-replay acceptance below.

## Bounded cache and memory

Default internal retained-checkpoint budget: **64 MiB and at most four checkpoints per filter instance**, whichever limit is reached first. No new public parameters. A checkpoint includes complete selected-plane L,C,Q arrays and metadata; logical state payload is 6*sizeof(float)=24 bytes per complex bin. Count actual allocated capacity/alignment/metadata toward the budget using checked arithmetic. Virtual S0 occupies no materialized cached arrays.

Publish a target S_n (n>=1) only after that request's output has completed successfully. Do not cache every replay step or retain a list of intermediate spectra. If a checkpoint exceeds the budget, skip caching; output still succeeds if required request workspace can be allocated. Cache capacity is an optimization, not admission based on clip length.

Use LRU eviction with synchronized lookup/access/publication. Duplicate publication of the same n may keep either identical state; do not mutate a live snapshot. A request lease keeps an evicted snapshot alive. Consequently the 64 MiB bound covers **cache-owned** storage, not all leased/request-private states: document peak as cache budget plus bounded state/scratch per simultaneously active request. Each request retains at most one starting checkpoint and one mutable full state, plus ordinary bounded transform/render scratch, immutable model leases and one acquired source frame. Reuse/move ownership when possible, without weakening isolation. Workspace retention is separately bounded by [execution](execution.md).

Neo-FFT/DS2-owned replay state and source retention are independent of clip length and replay distance for fixed geometry/concurrency. Host-held dependency pixel references require separately verified effective early release. Host history metadata, upstream caches and frames retained through properties are outside this storage guarantee and must be reported separately; no unconditional total-host-memory bound is claimed. Worst-case seeking still requires O(n) work from S0; finite checkpoints cannot guarantee constant replay latency for arbitrary requests. Do not hide this by capping replay distance, returning stale state or allocating clip-proportional memory. Publish diagnostics for replay start, steps, cache hits/bytes and peak retained source owners through internal tests/benchmark output, not video properties.

## Failure and concurrency

Requests never mutate a shared current Kalman state. Concurrent forward/backward/duplicate requests may duplicate work but must produce canonical output. A failure at an intermediate frame, FFT, model, output allocation or conversion discards private work and publishes no partial target checkpoint. Existing checkpoints remain valid and future retry starts from a valid snapshot. No poisoned ready flag or lost lease.

Do not hold cache/model/pool locks across host callbacks, frame acquisition, FFT or worker waits. Locks protect short ownership transitions only. Teardown waits for/releases active leases according to the host lifecycle, then destroys states, models, plans and pools without referencing a destroyed environment. Repeated creation/destruction, downstream cancellation and delivery errors are acceptance cases.
