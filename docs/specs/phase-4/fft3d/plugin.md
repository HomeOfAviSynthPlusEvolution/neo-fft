# FFT3D phase-4 interface

Specification: F3D-VS-004. Extends [phase-3 interface](../../phase-3/fft3d/plugin.md). No signature changes.

| Parameter | Domain/default | Meaning |
| --- | --- | --- |
| bt | int, default 3; -1..5 | 0 opens Kalman; -1/enhancement and 1..5 remain inherited |
| kratio | binary32, default 2; >=0 | Kalman per-component motion threshold ratio; raw domain always checked |
| l,t,r,b | int32, default 0; >=0 | Full-resolution pixel margins; [ROI rules](kernel-roi-fields.md) |
| interlaced | bool, default false | Pack/unpack ROI rows using the specified full-height permutation |
| opt,ncpu,mt,measure | Existing types/defaults | [Execution mapping](../execution.md); mt=true and historical opt aliases opened |

Resolve precedence before allocating derived tables:

1. Parse all raw types/domains, finite conversions and inherited static constraints.
2. Resolve selected planes and ROI/field geometry. These apply even to preview, bt=-1 and bt=0 frame 0; invalid static geometry fails creation.
3. Apply the phase-3 effective-preview predicate unchanged. Effective pshow overrides bt, noise strengths, Kalman state and enhancement. Request only n, preview inside the ROI, copy its exterior; no state/cache initialization or replay.
4. Without preview, bt=-1 keeps enhancement-only behavior with ROI/packing. bt=1..5 keeps the phase-2 neighborhood/fallback and phase-3 selected noise/enhancement model with ROI/packing applied to every consumed frame.
5. Without preview, bt=0 uses [Kalman](kernel-kalman.md) and [replay](../replay.md). No temporal neighborhood/fallback and no temporal noise multiplier. beta is validated but not consumed. kratio is consumed only here.

For bt=0, frame 0 is a bitwise copy of all source visible planes and properties, including selected planes. It performs no transform, noise sampling, enhancement or source-content finite scan. It does not seed state from source frame 0. Static configuration is still validated at creation. A one-frame clip therefore has a valid all-copy Kalman output.

For n>0, retain source properties from n; unselected planes and selected-plane pixels outside ROI are bitwise copies from n. Reconstruct/convert only selected ROI pixels. Intermediate replay reads only selected ROIs; NaN/Inf outside those regions does not fail evaluation. Source dependencies that the host cannot deliver may still fail the request.

## Active constants

Use phase-3 noise selection, with this explicit exception: **Kalman always consumes scaled sigma1 for the initial covariance**, including positive-pfactor sampled mode. Compute R0=sigma_eff^2/norm per selected plane and require it finite at creation. Raw sigma2..4 remain checked; sampled mode does not derive their unused powers. For pfactor=0, the inherited scaled-sigma comparison selects uniform versus analytic noise. Pattern tables and sample coordinates are active only for effective Kalman n>0 processing, with static selection-grid validation at creation for that mode.

kratio^2 must be finite when Kalman is the effective mode. For uniform/analytic models validate the resulting motion thresholds at creation; sampled-model thresholds are validated during model construction. Positive pfactor scales sampled power as in phase 3, before the Kalman floor. Do not calculate Wiener beta gains or degrid corrections inside the recurrence. Build enhancement constants only for enabled enhancement; enhancement's separate degrid path remains inherited.

Invalid margins, unaligned selected-chroma margins, empty/inadmissible selected ROI, odd selected ROI height with interlaced=true, invalid bt and invalid kratio fail creation. Allocation failure is a resource error, not a request to substitute bt=1, disable a plane or shorten history. No new clip-length cap is introduced; use checked wide index arithmetic and on-demand stages instead of a clip-sized dependency array.
