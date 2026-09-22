# FFT3D canonical Kalman operator

Specification: F3D-KALMAN-004. Applies only to effective bt=0, n>0. Windows, spatial R2C layout and normalized inverse/reconstruction follow phase 1; ROI/row mapping follows [ROI and fields](kernel-roi-fields.md).

## State and initialization

For each selected plane, each spatial block in canonical Y-then-X order, and each logical spatial half-spectrum bin, retain complex L (last estimate), C (covariance) and Q (process covariance). Real/imaginary C,Q are separate scalar quantities, not complex arithmetic. Padding is not a state bin. No temporal FFT is used.

Let `norm=1.0f/(float(bw)*float(bh))`, the phase-1 inverse spatial transform volume. Uniform noise power is `R0=(sigma_eff*sigma_eff)/norm`; do not replace transform volume with window energy. Use the inherited format scaling and preserve this binary32 operation order. Virtual state S0 is:

```
L.re = L.im = 0
C.re = C.im = Q.re = Q.im = R0
```

This remains true for analytic/sampled noise: do not initialize C,Q with the per-bin pattern. Source frame 0 is never transformed to create S0. Output 0 is pass-through and does not advance the state. S_i for i>=1 means the state **after** consuming source frame i. All selected planes use their own R0/geometry.

## One step

Gather the current source ROI into its reflected (optionally field-packed) cover, subtract the inherited sample midpoint, apply the normal analysis window and compute the spatial FFT X. No degrid subtraction is applied to X for the recurrence.

For a uniform model use R=R0. For an analytic or sampled pattern use R=max(P[k],1e-15f), where P is the immutable spatial power table; it repeats across blocks. Analytic P is inherited unchanged. Sampled P is the unscaled measured power `|X_sample-gG|^2*w` from phase 3: any positive pfactor selects sampling, but its magnitude does not multiply P for Kalman. Construct this unscaled table directly; do not recover it by dividing an already scaled table, which would introduce rounding, underflow or overflow dependence on pfactor. Do not multiply by a temporal length. Let K2=kratio*kratio and H=R*K2. All are binary32; require finite derived R,H. The motion predicate is:

```
dr = X.re - L.re
di = X.im - L.im
motion = dr*dr > H || di*di > H
```

It is not squared complex magnitude. Equality takes the smoothing branch. Either component triggers a reset of **both** components. A squared difference overflowing to +infinity compares greater than a finite H and therefore resets; it need not become an error. Non-finite source spectra or persistent state are errors.

If the uniform R is zero, first take the explicit identity step `L=X; C=Q=0` for both components. This also handles X=L=0 without division by zero. Pattern models do not take that branch because of the floor.

Otherwise, if motion is true:

```
L = X
C.re = C.im = Q.re = Q.im = R
```

Otherwise, independently for real and imaginary components, using the old C,Q,L:

```
sum = C + Q
gain = sum / (sum + R)
Q_new = (gain * gain) * R
C_new = (1.0f - gain) * sum
L_new = (gain * X) + ((1.0f - gain) * L)
```

Preserve the shown binary32 operation order. The scalar oracle uses separate operations, without contraction/reassociation. Require finite active intermediates and new state; source-dependent overflow fails the request, not a partial checkpoint. In SIMD, reset/zero-noise lanes must not turn unused speculative divides into errors or write NaNs. Odd bin counts and tails access only logical elements.

## Output and state separation

At the requested n copy L_n to a private output spectrum, apply inherited sharpen/dehalo (including their specified degrid handling), then normalized inverse FFT and ordinary synthesis/crop/conversion. Neither enhanced spectrum, reconstructed pixels, clipping nor quantization feeds back into L,C,Q. Replay-only intermediate frames need no enhancement, inverse FFT or output allocation. Final output and state must match a canonical sequential run with enhancement enabled only at its output stage.

Example oracle: R0=R=2, kratio=2, X1=(1,0). No motion; initial sum=4, gain=4/6, so L1.re=2/3, Q1.re=(4/6)^2*2 and C1.re=(1-4/6)*4 with binary32 rounding. With X1=(3,0), 9>8 resets both components to L1=(3,0), C1=Q1=(2,2). These distinguish zero initialization, component-wise reset and a mistaken seed from frame 0.

An implementation may vectorize independent bins/blocks/planes. It may not parallelize dependent time steps, merge covariance components, feed enhancement back or replace early history with an arbitrary fresh state. Same-build checkpoint restore must retain exact binary32 state, not compressed, quantized or reconstructed state.
