# DFTTest temporal overlap-add operator

Specification: DFT-TOLA-005. Applies to tmode=1. Let B=sbsize, T=tbsize, O=tosize, D=T-O, V=T*B*B, K=floor(B/2)+1. Reuse phase-1 spatial origins/padding and the phase-2 [3D transform layout](../../phase-2/kernel-fft.md): real [T,B,B], compressed [T,B,K], x contiguous. Even T is a full complex temporal axis, not a second compressed axis.

## Filtering window

Compute raw binary64 temporal samples `r[z]=raw_window(z,T,twin,tbeta)` using the inherited half-sample equations. Retain the raw array until all denominators have been computed. For each z:

```
e[z] = 0.0                         // binary64
for h=z; h>=0; h-=D: e[z] += r[h]*r[h]
for h=z+D; h<T; h+=D: e[z] += r[h]*r[h]
tw[z] = r[z] / sqrt(e[z])
```

Require e[z] finite and >0 and tw[z] finite. Preserve this summation order rather than assuming all members of a residue class round identically. Do not use already-normalized values in later denominators. Even O=0 normalizes: `tw[z]=r[z]/abs(r[z])` when nonzero, not the raw window. A zero phase energy is a creation error, not epsilon substitution. Keep the sign of negative raw window samples.

Compute sw by the phase-1 spatial rules (normalize only for smode=1), then:

```
h[z,y,x] = float(tw[z]*sw[y]*sw[x]*(1/sqrt(V)))
E = sum_float_in_z_y_x_order(h*h)
wscale = 1.0f/E
```

Use binary64 window products followed by one storage cast. Require finite positive E/wscale. Energy comes from the actual stored float window, not its raw temporal version or only the target slice. There is no extra factor for number of overlapping blocks. The normalized h is shared by all temporal starts for the same configuration/plane geometry.

## Spectrum, mean and models

At each temporal start s and each spatial origin, gather the T logical slots and multiply inherited sample amplitudes by h. Run a full 3D forward FFT. Phase-2 zmean applies once per entire volume using `G=FFT(255*h)` and one DC ratio; do not remove means per frame/slice. For zmean=true require finite G and nonzero real DC. Restore the window-shaped mean before inverse.

All five filters keep phase-2 PSD and gain definitions: power is re^2+im^2 without Hermitian doubling or T scaling. Recompute their consumed constants using this window's wscale. Scalar sigma for ftype<2 is divided by wscale; gain-type sigma is not. pmin/pmax retain their inherited scaling. No second sigma squaring or overlap-count multiplier.

Phase-3 curves use the same [T,B,K] frequency coordinates, including even T's Nyquist plane. Active nlocation still replaces the primary scalar/curve table. Its h2/G2 use **raw temporal and spatial windows with both overlap normalizations disabled**, independent of output tmode/smode/O. The existing sampled-power calibration `A=Psum*((1/M)*(wscale2/wscale)*alpha)` now uses this tmode=1 filtering wscale; no extra division, sqrt or O/D factor. Sampling slots remain fn..fn+T-1, not lattice-aligned output blocks. Validate inactive raw arrays but do not calculate overridden derived models.

## Ordered reconstruction

Initialize one request-private output-plane accumulator to binary32 zero. Process covering temporal starts in ascending order. Within each start process spatial origins in canonical Y-then-X order. After mean restoration let U be the inverse normalized by 1/V exactly once. At target slice z=n-s compute:

```
contribution[y,x] = (U[z,y,x]*float(V))*h[z,y,x]
```

For smode=1 add all local B*B samples at that spatial origin. For smode=0 add only the local spatial center (B/2,B/2) into its output position. Temporal contributions use `+=` in **both** spatial modes; assigning the smode=0 center as in tmode=0 would discard earlier time blocks.

Do not pre-sum each temporal block into a separate rounded image and then add images unless it reproduces the canonical contribution order exactly. Parallel transforms are allowed; commit their contributions in that order. Once every covering start has contributed, crop and perform inherited conversion/dither once. No per-block clipping, quantization or diffusion; no scaling by number of contributors. Only n is published, even if an optimized inverse computes other slices.

For an identity spectrum, ideal temporal gain is `sum_s tw[n-s]^2=1` at all n, including the endpoints. smode=1 adds its inherited spatial squared-overlap gain; smode=0 retains its spatial center factor `sw[B/2]^4`. The stored-float oracle is the ordered sum of `V*h[z,local_y,local_x]^2` over contributing origins. Do not force exact unity by dividing the output by a measured weight sum.

T=1,tmode=1 still normalizes tw[0]. Unlike tmode=0, this removes any non-unit raw singleton temporal magnitude. Use its exact table/reconstruction contract rather than a blanket bypass or cross-mode bit-equality assertion.

## Small oracles

- B=1,smode=0,rectangular windows,T=4,O=2: tw=1/sqrt(2), h=1/sqrt(8), E=1/2,wscale=2. Two covering blocks each contribute one-half of the identity sample. With zmean=true,ftype=2,sigma=0, each contributes one-half of its four-slot mean. At n=2 the result is `(mean(q[0:4])+mean(q[2:6]))/2`, not the centered-mode mean.
- B=1,rectangular windows,T=5,O=2: phase classes 0 and 1 have two members, class 2 one. tw=[1/sqrt(2),1/sqrt(2),1,1/sqrt(2),1/sqrt(2)], ideal E=3/5,wscale=5/3. n=2 has one contribution of weight 1; n=3 has two of weight 1/2.
- O=0: one covering block contributes to each output; nonzero raw temporal samples normalize to signs and temporal identity gain is one. Non-identity filtering still depends on window signs and the full time block.

Examples allow prescribed floating rounding. They distinguish block phase, phase-dependent overlap counts, energy calibration and use of the target slice.
