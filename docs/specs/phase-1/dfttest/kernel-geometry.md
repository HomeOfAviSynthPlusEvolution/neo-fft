# DFTTest spatial block geometry

Specification: DFT-GEO-001. Evidence: static extraction from the fixed geometry, padding and CPU batch pipeline; safety admission is a new decision.

## Inputs and outputs

Inputs: actual selected-plane W,H, B=sbsize>=1, smode in {0,1}, and sosize. Temporal extent is 1. Outputs: effective overlap O, step S, padded dimensions P,Q, source offset dx,dy and ordered block origins. Follow [MEM-001](../kernel-plane-geometry.md) for checked integer arithmetic, views and aliasing.

For smode=0, require odd B, and set O=0 regardless of the supplied sosize. For smode=1 set O=sosize, require 0<=O<B, and if O>B/2 require B divisible by B-O. This restriction applies even though other overlaps could be mathematically normalized.

## Center-sample mode: smode=0

Let c=floor(B/2):

```text
P=W+2*c; Q=H+2*c
dx=dy=c; S=1
origins: y=0..H-1, x=0..W-1, Y outer then X inner
```

Each block is BxB. Emit its local sample (c,c) to padded accumulator position (y+c,x+c), then crop at (dx,dy). Equivalently it directly produces visible output (y,x). Every output has exactly one selected center, not a sum over all blocks containing it.

## Overlap-add mode: smode=1

```text
S=B-O
A=2*max(S,O)
P=B*ceil(W/B)+A
Q=B*ceil(H/B)+A
dx=floor((P-W)/2)
dy=floor((Q-H)/2)
Nx=floor((P-B)/S)+1
Ny=floor((Q-B)/S)+1
```

Origins are (i*S,j*S), i=0..Nx-1 and j=0..Ny-1, in Y-then-X order. Equivalently the reference's Y loop is `y < floor((Q-O)/S)*S` with step S. Only complete blocks are transformed. A trailing cover strip beyond the last block can exist; the visible crop must have full phase coverage and never include an unwritten tail.

Each block contributes every local sample to its corresponding padded accumulator position. After all contributions, crop `[dy:dy+H,dx:dx+W]`. P and Q are cover dimensions, not rounded output dimensions.

## Reflection and storage

Both modes require dx<=W-1, P-dx-W<=W-1, dy<=H-1 and Q-dy-H<=H-1. Fill the padded source using the endpoint-excluding reflection in MEM-001. This admission applies per selected plane, so a luma-valid configuration may still be invalid on subsampled chroma.

One block has N=B^2 real values and B*(floor(B/2)+1) stored complex bins. Row pitch/alignment and reference scratch padding are not algorithm semantics. Unselected planes require only copying and do not require valid block padding; configuration enums and deferred features are still validated globally.

## Examples and failures

- W=5,H=4,B=3,smode=0 yields P=7,Q=6,dx=dy=1 and 20 block centers. A single visible impulse affects neighboring blocks' contents, but each output is emitted only by its own centered block.
- W=17,H=13,B=8,O=6,smode=1 gives S=2,A=12,P=36,Q=28,dx=9,dy=7,Nx=15,Ny=11. Right/bottom pads are 10,8 and satisfy one-reflection admission.
- B=8,O=5 is rejected because 8 is not divisible by 3. B=8,O=4 and B=8,O=0 are admitted overlap values subject to plane geometry.
- B=1,smode=0,1x1 source has no padding and is valid. B=1,smode=1 adds a one-pixel pad on each side and fails on a 1x1 source. Do not turn the failure into repeated reflection implicitly.

Reject invalid overlaps/modes, unsafe reflection, uncovered visible pixels and overflow at creation. Test first/last block and crop coordinates, non-divisible dimensions, odd/even B, both modes, heavy overlap and subsampled planes independently of FFT results.
