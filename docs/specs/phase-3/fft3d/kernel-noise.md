# FFT3D: frequency-dependent noise and preview

Specification: F3D-NOISE-003. Uses phase-1 spatial FFT/windows/grid G and phase-2 temporal filtering. W=bw, H=bh, K=floor(W/2)+1, norm=1/(WH). Valid bins are y=0..H-1, x=0..K-1. Tables have H*K logical floats; padding contributes to no reduction. Table arithmetic and row-major reductions are binary32; reject non-finite intermediates.

## Analytic sigma profile

Scale all sigmas by 2^(bits-8) for integer input, 1/255 for float, as in phase 1:

```text
fy = (H - 2*abs(y-floor(H/2))) / float(H)
fx = x / float(K)
f = sqrt((fx*fx + fy*fy)*0.5)
a = sqrt(0.5)/4; b = sqrt(0.5)/2
s = sigma4 + (sigma3-sigma4)*f/a                 if f<a
    sigma3 + (sigma2-sigma3)*(f-a)/(b-a)        if a<=f<b
    sigma + (sigma2-sigma)*(1-f)/(1-b)          otherwise
P[y,x] = s*s/norm
```

Do not substitute K-1 or periodic min-distance for these coordinates. For odd H, this reference fy is not symmetric about DC. Equality at knots uses the following branch; the function is continuous. This two-dimensional table is shared across temporal frequencies after the scaling below.

## Sampled profile

Use the normal cover, midpoint subtraction, analysis windows and forward FFT on the sample frame. For each block spectrum X compute g=degrid*Re(X[0])/Re(G[0]), or exactly 0 when disabled; keep g fixed for the whole block. R=X-gG:

```text
fy = 2*min(y,H-y)/float(H); fx = 2*x/float(W)
w[y,x] = (fy*fy+fx*fx)/(fy*fy+fx*fx+pcutoff*pcutoff)
score(X) = sum_yx (Re(R)^2+Im(R)^2)*w[y,x]
```

The sum has no conjugate multiplicity weights and no epsilon; w[0,0]=0. Manual (px,py) indexes the grid, at cover origin (px*step_x,py*step_y). Automatic (0,0) searches by=2..Ny-3, bx=2..Nx-3 in Y-then-X order, taking the first minimum with strict `<`. Initialize from the first candidate, not a numerical sentinel. A non-finite candidate fails the frame, rather than being skipped.

For the selected block:

```text
P = pfactor * |X-gG|^2 * w
```

No pfactor square, sigma substitution, sqrt, extra FFT-volume scale or window-energy correction is allowed. This is already measured spatial FFT power. Optional internal diagnostic psigma of the unscaled profile is sqrt(sum(|R|^2*w)/(sum(w)*H*K)); it does not filter or round the table.

## Consumption

For btcur=T, including fallback T=1, per-temporal-bin noise is T*P[y,x]. Uniform mode retains T*sigma_eff^2/norm. Replace only the noise term in [phase-2 Wiener](../../phase-2/fft3d/kernel-temporal.md); retain epsilon, beta floor and temporal-DC degrid/restoration. Do not square P again or divide it by T. Validate finite multipliers for every admitted T before model publication.

Store base P, not a mutable pattern3d initialized using the first request's T. Boundary-first evaluation must not freeze T=1 calibration for later interior frames.

## pshow image

Select a block from frame n using normal immutable analysis windows. Gather n again using rectangular analysis windows (all weights 1), retain only the selected block's raw spatial spectrum, and zero all other block spectra. Do not degrid, Wiener-filter or enhance the preview spectrum. Invert with normal spatial normalization and accumulate/crop using rectangular synthesis weights too. Restore the usual plane midpoint and perform ordinary conversion/clipping. Outside the selected block's cropped footprint the accumulator is zero; integer chroma returns its midpoint there. Only the chosen block contributes in overlap regions; do not introduce an overlap divisor.

Both gathers use normal plane IsChroma semantics. No persistent window changes, pframe fetch, temporal fetch or persistent sampled model is required. No labels are drawn (reference DrawString is disabled). Repeated/parallel/descending requests must match isolated requests.
