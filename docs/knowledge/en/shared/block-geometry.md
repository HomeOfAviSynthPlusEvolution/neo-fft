# Blocks, overlap and boundaries

[Contents](../README.md)

The following describes one axis; apply it independently vertically. W is the processed plane length, B block length, O overlap and S step. Chroma uses its own dimensions, but block sizes are not automatically reduced by subsampling.

## FFT3D

Negative O becomes floor(B/3); otherwise require `0≤O≤floor(B/2)`, B≥2:

```text
S = B-O
N = ceil((W-O)/S)+2
P = N*S+O
offset = S
block_origin(i) = i*S, 0<=i<N
```

W is the selected plane ROI length, at least B, with legal reflection as well. W=64,B=8,O=2 gives S=6,N=13,P=80, left pad 6 and right pad 10. The visible interval is `[6,70)`. ROI limits filtering, not output dimensions. When chroma is selected, luma-unit margins must be divisible by the relevant subsampling factor.

## DFTTest

Center mode smode=0 requires odd B. With c=floor(B/2), P=W+2c, offset=c and S=1; there are W blocks, each emitting only its local center. W=5,B=3 gives P=7 and five centers corresponding to five output samples.

Spatial overlap mode smode=1 uses:

```text
S = B-O
A = 2*max(S,O)
P = B*ceil(W/B)+A
offset = floor((P-W)/2)
N = floor((P-B)/S)+1
```

Require `0≤O<B`; when O>floor(B/2), B must be divisible by S. W=17,B=8,O=6 gives S=2,P=36,offset=9,N=15. An unused trailing cover strip may exist, but the visible crop must be fully covered. This is not FFT3D's grid formula.

## One reflection

Subtract offset from a padded coordinate to obtain j. Read:

```text
j < 0       -> -j
0 <= j < W  -> j
j >= W      -> 2*W-2-j
```

Padding `[a,b,c,d]` by two samples on both sides yields `[c,b,a,b,c,d,c,b]`, without repeating endpoints. Neither pad may exceed W−1. Invalid geometry fails creation instead of repeatedly applying reflection. Small chroma planes can therefore fail before luma does.

FFT3D field packing orders a height-8 ROI as `[0,2,4,6,7,5,3,1]`, applies the geometry to the entire packed image and inverts the permutation afterward. It does not process two separate half-height images. DFTTest temporal endpoints clamp rather than reflect; see [temporal overlap](../dfttest/temporal-ola.md).

Integer geometry is checked for overflow before allocation. Increasing overlap increases block count and spectral work without changing output resolution. See [FFT3D computation](../fft3d.md) and [DFTTest API](../../../api/en/dfttest.md).
