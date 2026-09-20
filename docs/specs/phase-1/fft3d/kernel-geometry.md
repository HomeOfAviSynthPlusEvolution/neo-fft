# FFT3D spatial block geometry

Specification: F3D-GEO-001. Evidence: static extraction from the fixed FFT3D engine/frame-cover code; explicit safety admission is a new decision.

## Inputs and outputs

Inputs: actual selected-plane width W and height H, block Bx=bw, By=bh, overlaps ow,oh. Phase 1 requires Bx,By >= 2, W>=Bx and H>=By. Output: effective overlaps, block counts, cover dimensions, source offset, block origins and crop. All dimensions/offsets satisfy [MEM-001](../kernel-plane-geometry.md).

An overlap parameter less than zero means floor(B/3), including values below -1. Otherwise require `0 <= O <= floor(B/2)`. Let Sx=Bx-Ox and Sy=By-Oy. Both are positive. Dimensions and overlaps are per-plane sample counts, not luma-scaled values.

## Grid and extension

```text
Nx = ceil((W-Ox)/Sx) + 2
Ny = ceil((H-Oy)/Sy) + 2
P  = Nx*Sx + Ox
Q  = Ny*Sy + Oy
dx = Sx
dy = Sy
```

Block (i,j), 0<=i<Nx and 0<=j<Ny, starts at (i*Sx,j*Sy) and contains Bx*By samples. The last block ends exactly at the cover edge. Source pixels occupy `[dx,dx+W)` by `[dy,dy+H)`.

Before constructing the cover, require each axis to meet the one-reflection admission in MEM-001: dx and P-dx-W are at most W-1; similarly dy and Q-dy-H at most H-1. Fill `cover[y,x]=source[r(x-dx,W),r(y-dy,H)]`. Reflection excludes endpoint duplication and is progressive; interlaced rows are not supported.

The output is the reconstructed cover crop `[dy:dy+H, dx:dx+W]`. There is no visible output padding, size rounding, block-aligned crop or changed chroma size. Unselected planes do not need covers.

## Storage and examples

Each real block has Bx*By samples; its spectrum has By*(floor(Bx/2)+1) complex bins. Physical alignment/pitch and batch capacity are private; old pitch rounding to 16 elements has no output semantics. Process blocks in row-major origin order when accumulating. FFT input/output storage is disjoint.

- W=64,H=48,Bx=8,By=6,Ox=Oy=2 gives steps 6,4, counts 13,14, cover 80x58, source offset (6,4), and right/bottom pads 10,6. There are 182 transforms of shape [6,8].
- B=8 and omitted O gives O=2, not 8/4 as a separate rule. B=32 gives O=10 and S=22.
- W=8,B=8,O=0 gives Nx=3,P=24,dx=8; dx>W-1, so phase 1 rejects this geometry. Do not run a reference with unsafe padding to create a golden output.
- For YUV420, luma 128x96 and block 8x6 use the same 8x6 block on the actual 64x48 chroma planes, with geometry derived separately.

## Errors and acceptance

Invalid blocks, overlaps, cover padding, arithmetic overflow and unrepresentable storage are creation errors naming the affected plane. Verify odd/even and non-square blocks, zero and half overlap, negative-overlap defaulting, non-divisible image dimensions, subsampled plane boundaries and the first accepted/rejected reflection sizes. A tiny-plane rejection is an intentional scope restriction, not a successful comparison of old undefined behavior.
