# DFTTest absolute temporal-block geometry

Specification: DFT-TGRID-005. Inputs are admitted T=tbsize, O=tosize, N=clip length and output index 0<=n<N, with tmode=1. Let D=T-O and A=-max(D,O). All index arithmetic below uses checked signed wide integers; floor division means mathematical floor, not C++ truncation toward zero.

## Fixed lattice and covering blocks

Temporal block starts are `s_j=A+j*D`, j>=0. Only blocks with `s_j<=n<s_j+T` contribute to output n:

```
j_first = max(0, floor_div(n-T-A,D) + 1)
j_last  = floor_div(n-A,D)
for j = j_first .. j_last, ascending:
    s = A + j*D
    z_target = n-s                 // 0 <= z_target < T
```

There is at least one covering block, and at most ceil(T/D). Generate starts directly in O(ceil(T/D)) work; do not copy the original's loops from the negative anchor to n. Do not reset the lattice at the first requested frame, at a cache miss, a worker boundary or a scene cut.

The negative anchor matters, especially at low overlap: T=5,O=2 gives A=-3,D=3, while T=5,O=4 gives A=-4,D=1. The phase is not generally equivalent to starting at zero. No real output before 0 or after N-1 is created, but blocks can extend beyond either end.

For each covering block preserve T logical source slots:

```
logical_time(z) = s+z, z=0..T-1
source_index(z) = clamp(s+z,0,N-1)
```

Clamping applies to each slot after fixing the block start, not to s itself. Duplicate source indices retain distinct positions/windows in the transform. In particular, different starts near a boundary remain different blocks even if some/all real indices coincide. Never key a cached block by a clamped start or only by its deduplicated source set.

The target slot always maps to source n. Every consumed slot is a whole source frame before inherited per-plane spatial padding; no temporal reflection, wrapping, field extraction or shortened transform. T<=N is retained even though clamping would make some longer transforms mathematically possible.

## Examples

| T,O,n | Covering starts s | Target positions n-s |
| --- | --- | --- |
| 4,2,0 | -2,0 | 2,0 |
| 4,2,3 | 0,2 | 3,1 |
| 5,2,0 | -3,0 | 3,0 |
| 5,2,2 | 0 | 2 |
| 5,2,3 | 0,3 | 3,0 |
| 4,0,3 | 0 | 3 |
| 4,0,4 | 4 | 0 |
| 6,4,0 | -4,-2,0 | 4,2,0 |

With N=5,T=4,O=2,n=4, starts are 2,4; the logical slot maps are [2,3,4,4] and [4,4,4,4], with target positions 2,0. Use their different window slots even though both include repeated last-frame samples.

All source indices in contributing blocks lie in the clamped range [n-(T-1),n+(T-1)]. Their deduplicated union therefore contains at most min(N,2*T-1) real frames (at most 29 here). This is a finite local dependency set, not input-history recurrence. Additional active noise-sample dependencies are specified separately in execution.md.

The fixed original's cold ring initialization may transform older blocks not covering n. Omit those and their source fetches. Reference comparisons concern their contribution to output n, not matching internal request counts or failures on unused samples.
