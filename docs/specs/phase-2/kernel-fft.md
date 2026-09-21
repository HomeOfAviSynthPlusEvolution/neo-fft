# Phase-2 multidimensional real FFT adapter

Specification: FFT-TIME-001. Extend [FFT-001](../phase-1/kernel-fft.md), retaining binary32 real storage, live complex<float> objects, read-only input, disjoint input/output, checked descriptors, private destructive-backend scratch and normalized inverse.

## Shape and batch

DFTTest real shape is [T,H,W], compressed complex shape [T,H,K], K=floor(W/2)+1. Logical order is z,y,x with x fastest (frequency order kt,ky,kx). A batch contains independent **whole volumes**, never independent temporal slices. Real element count is V=T*H*W, complex count T*H*K; twice the latter is a float-component count, not a complex count.

Descriptors specify element strides for each axis, outer batch distance, capacity C and active count a. Dense real strides are [H*W,W,1], complex [H*K,K,1]. Non-dense positive, non-overlapping layouts require checked final addresses and sufficient storage. Shape products, float/complex byte conversions, stride offsets and batch capacities must not overflow. No implicit padded real columns or vector alignment beyond natural alignment are required of callers.

Only active volumes and logical elements may be touched. a=0 is a no-access no-op; padding and inactive slots remain untouched. T=1 may use the 2D backend plan provided the logical definition is unchanged.

## Transform and half-spectrum

```text
X[kt,ky,kx] = sum(z=0..T-1,y=0..H-1,x=0..W-1)
  a[z,y,x]*exp(-2*pi*i*(kt*z/T + ky*y/H + kx*x/W))

a[z,y,x] = (1/V)*sum(all full kt,ky,kx)
  X[kt,ky,kx]*exp(+2*pi*i*(kt*z/T + ky*y/H + kx*x/W))
```

Forward scale is 1; inverse scale is 1/V, excluding batch count. Missing coefficients are reconstructed by simultaneous conjugation of **all three axes**:

```text
X[kt,ky,kx] = conjugate(X[(-kt) mod T,(-ky) mod H,(-kx) mod W])
```

kx=0 and the even-W Nyquist column are not generally real. Only coordinates self-conjugate on every axis have zero imaginary part. Odd W has no Nyquist column; odd T has no temporal Nyquist plane. A 2D spatial conjugation rule applied independently to every kt slice is wrong.

Parseval check: sum(a*a)=sum(weight[kx]*abs(X)^2)/V over stored bins, with weight=1 for kx=0 and even-W Nyquist, 2 otherwise. These weights verify total energy only; neither filter doubles per-bin PSD or thresholds. Symmetric real gains and valid mean templates preserve Hermitian compatibility. Arbitrary incompatible spectral test input is outside C2R's domain and must fail validation rather than define new output behavior.

## Relationship to FFT3D

FFT3D retains the 2D real adapter for each frame and performs a small **complex** temporal DFT per spatial bin as defined in [F3D-TIME-001](fft3d/kernel-temporal.md). It keeps all T temporal coefficients and selectively reconstructs cur. Its 1/T temporal inverse and 1/(bw*bh) spatial inverse are separate; do not reuse DFTTest synthesis scaling without adapting the algorithm.

## Required independent checks

Use binary64 direct DFT to check forward sign, axis ordering and normalization, not merely a forward/inverse round trip. Include non-square adapter shapes [3,2,5] and [3,3,4], spatial/temporal impulses, constants, sinusoidal input, strided volumes, active counts 0,1,C-1,C, guard padding and T=1 equivalence. A unit impulse at z=1,y=x=0 has X[kt,ky,kx]=exp(-2*pi*i*kt/T). A constant q has DC=V*q and normalized inverse q. No batch-count factor may enter either value.

No global fast-math/implicit contraction or hidden FFT worker pool is introduced. Backend revision, threading, plan ownership, failure categories and concurrency rules remain those of phase 1.
