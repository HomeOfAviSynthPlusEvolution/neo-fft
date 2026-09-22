# DFTTest: phase-3 VS interface

Specification: DFT-VS-003. Extends [phase-2 interface](../../phase-2/dfttest/plugin.md). Temporal/window/mode/format/execution restrictions remain unless opened here.

| Parameter | Default | New domain/effect |
| --- | --- | --- |
| slocation | empty float array | Frequency/value pairs for shared curve; nonempty overrides axis arrays |
| ssx, ssy, sst | empty float arrays | X, Y, time curves when slocation is empty |
| ssystem | 0 | 0 separable, 1 radial; exact legacy selection below |
| nlocation | empty int array | Up to 500 (start_frame, plane, y, x) tuples |
| alpha | 5 for ftype=0, otherwise 7 | Finite >0; estimated power multiplier |

All float array elements must be finite and representable in binary32. Each nonempty curve must have even length, include endpoint pairs at 0 and 1, have positions in [0,1], nonnegative values, and unique positions after conversion. Arbitrary input order is allowed; sort by position. Validate ignored axis arrays too. Do not insert endpoints, average duplicates or clamp positions.

nlocation length must be divisible by 4, with <=500 tuples. For each require 0<=start_frame<=N-T, an existing plane and 0<=x<=plane_width-S, 0<=y<=plane_height-S, where S=sbsize. Coordinates refer to actual plane pixels (RGB planes are not subsampled). No padding/reflection, endpoint clamping, centering or coordinate scaling applies to sample rectangles. Validate even with planes=[] or ftype>=2. Repeated tuples are legal and retain multiplicity in the average.

For ftype=0/1 and nonempty nlocation, the estimated table **replaces** the scalar/curve primary sigma table, rather than adding to or multiplying it. For ftype=2/3/4, nlocation/alpha are validated but inactive and cause no dependencies. sigma2/pmin/pmax remain scalar with existing defaults/scaling. Curves may provide primary sigma for all five ftypes; never square values just because they are called sigma.

A tuple may name an output-unselected plane; that patch is explicitly consumed to form the shared model. All selected output planes use the same estimated table, not separate plane estimates. With planes=[], validate configuration then copy n only without estimation. See [profile contract](kernel-profile.md) for exact arithmetic and [execution](../execution.md) for failures/publication.
