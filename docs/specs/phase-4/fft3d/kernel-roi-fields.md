# FFT3D ROI and field packing

Specification: F3D-ROI-004. Applies to all effective FFT3D modes except frame-0 Kalman output, which is a complete copy after static validation.

## Coordinates and admission

l,t,r,b are full-resolution margins. Require checked l+r<source_width and t+b<source_height. For a selected subsampled YUV chroma plane with horizontal/vertical subsampling factors ax=2^subW, ay=2^subH, require l,r divisible by ax and t,b divisible by ay. Divide those margins exactly; RGB, Y and GRAY use factors 1. Unselected chroma imposes no crop-alignment constraint because it is copied. Existing source-format/dimension validation still applies.

For each selected plane use its actual width Wp/height Hp and derived margins lp,tp,rp,bp:

```
W = Wp - lp - rp
H = Hp - tp - bp
source_roi(x,y) = source_plane(lp+x, tp+y)
```

Use W,H everywhere in phase-1 grid/cover admission, including W>=bw, H>=bh and one-reflection bounds. Never derive the block count from a differently rounded full-resolution crop. No padding or block may read original pixels outside the ROI; outer cover samples reflect the ROI itself with the inherited endpoint-excluding reflection. Validate all count/offset/byte products with checked arithmetic before narrowing/allocation.

Copy selected full visible planes from source n before writing the filtered ROI, or use an equivalent partitioned copy with identical visible bytes. Unselected planes are complete copies. No read/write outside a visible row through pitch padding. Output size/format/frame count/rate and properties remain from source n.

## Interlaced mapping

With interlaced=false, the working image is the ordinary ROI. With interlaced=true require **H even for every selected plane**, then pack an H-row working image before reflection/windowing:

```
q = H/2
source_row(v) = 2*v                  for 0 <= v < q
source_row(v) = 2*(H-1-v) + 1        for q <= v < H
packed(x,v) = source_roi(x,source_row(v))
```

For H=8, rows are `[0,2,4,6,7,5,3,1]`. After synthesis, unpack with the inverse permutation and write `(lp+x,tp+source_row(v))`. Margins may place tp on an odd absolute row; even/odd in this mapping is **ROI-relative**, matching the reference. No additional top-margin parity restriction is imposed. The selected ROI height, rather than unselected planes or full-frame height, determines admission.

Grid/cover/block height stays H. Do not halve bh, create two field grids or reflect each half separately. Blocks can cross the middle seam. Apply outer reflection to the already packed image. The permutation is applied independently and identically to every temporal/sample/replay frame. Frame indices remain full-frame indices; do not double temporal length or infer TFF/BFF from properties. Preserve `_FieldBased` and all inherited properties; this is neither deinterlacing nor a field-order conversion.

## Interactions

- Temporal FFT3D uses the same cropped/packed geometry for every logical slot; boundary fallback rules are unchanged.
- Pattern pframe and preview n use cropped/packed cover blocks. px,py remain **block-grid** coordinates, not original pixel coordinates. Automatic candidate admission uses the ROI-derived grid.
- pshow writes its preview only inside the ROI, using the inherited selection/preview windows and inverse row mapping. Exterior pixels are copied from n.
- Kalman retains one state for each bin of the full packed grid; packing does not create two recurrence timelines. Changing geometry requires a new filter instance and cannot reuse checkpoints.
- Integer YUV chroma midpoint subtraction and output restoration depend on plane identity/format, not crop origin, packed position or preview mode. Float handling is inherited without an invented midpoint.

Negative margins are rejected, not clamped as in the legacy implementation. Misaligned subsampled margins and odd packed heights are rejected rather than reproducing reference size disagreement or out-of-bounds access.
