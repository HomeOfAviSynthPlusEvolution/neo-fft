# FFT3D noise sources and enhancement order

[Contents](../README.md)

The [FFT3D overview](../fft3d.md) defines spatial/temporal Wiener filtering. This article explains its noise powers and why enhancement order matters. See [API domains](../../../api/en/fft3d.md). W=bw,H=bh; formulas traverse only logical half-spectrum bins.

## Three noise sources

At creation, pfactor>0 selects sampling; otherwise unequal effective sigmas select an analytic profile; otherwise use uniform noise. Scale sigmas by `2^(bits−8)` for integers or `1/255` for float. Uniform spatial power is `P=W*H*sigma_eff²`, using block area rather than window energy.

The analytic profile uses already-scaled sigmas:

```text
norm = 1/(W*H)
fy = (H - 2*abs(y-floor(H/2)))/H
fx = x/(floor(W/2)+1)
f = sqrt((fx*fx+fy*fy)/2)
a = sqrt(0.5)/4; b = sqrt(0.5)/2
s = sigma4 + (sigma3-sigma4)*f/a             if f<a
    sigma3 + (sigma2-sigma3)*(f-a)/(b-a)    if a<=f<b
    sigma + (sigma2-sigma)*(1-f)/(1-b)      otherwise
P = s*s/norm
```


Do not replace the fx denominator with floor(W/2) or fy with the usual symmetric frequency mapping, especially for odd H. Ordinary T-frame Wiener uses T*P. Stored base P must not depend on whether the first request falls back to one frame.

## Sampling, scoring and preview

Use the clamped pframe with ordinary reflection/analysis windows to obtain X. Remove the window-shaped grid using `g=degrid*Re(X[0])/Re(G[0])`; degrid=0 directly uses R=X. Define `fy=2min(y,H−y)/H`, `fx=2x/W`, `w=(fy²+fx²)/(fy²+fx²+pcutoff²)`. Score is Σ|R|²w without half-spectrum multiplicity.

Manual px/py are block indices. Both zero searches bx=2..Nx−3,by=2..Ny−3 in Y/X order, retaining the first minimum. Ordinary sampled power is `P=pfactor*|R|²*w`. For |R|²=100,w=0.4,pfactor=2, P=80 and T=3 uses 240. Kalman is different: positive pfactor selects sampling but does not scale the measured power, using 40 here.

When sampling or an unequal-sigma model is active, pshow overrides other paths. It selects a block from current n with normal windows, then reconstructs that block alone with rectangular analysis/synthesis windows. No denoising, enhancement, pframe fetch or sampled-model publication occurs. Other block contributions inside the ROI are zero; integer chroma restores its midpoint. This is not a denoised noise-map image.

## Frequency weights and enhancement gains

Recompute g from the spatial spectrum Z being enhanced and set R=Z−gG. smin/smax use sigma's format factor; ht **does not**:

```text
dy = y if y<floor(H/2) else H-y
d2 = dy*dy*svr*svr/floor(H/2)^2 + x*x/floor(W/2)^2
Ws = 1-exp(-d2/(2*scutoff*scutoff))
rawWh = exp(-0.7*d2*hr*hr)-exp(-d2*hr*hr)
Wh = rawWh/max(rawWh)
A = smin_eff^2/norm; B = smax_eff^2/norm; C = ht^2/norm
q = abs(R)^2+1e-15
S(q) = 1 + sharpen*Ws*sqrt(q*B/((q+A)*(q+B)))
D(q) = (q+C)/((q+C)+dehalo*Wh*q)
output = R*S(q)*D(q)+gG
```


Zero strengths bypass their gains exactly; dehalo=0 skips unused Wh normalization. Active Wh needs a positive finite maximum. B=0 makes sharpen gain 1. For H=5, dy is `[0,1,3,2,1]`, not symmetric minimum distance.

Both gains use the same residual power. For q=100,A=16,B=400,C=2500,Ws=Wh=1,sharpen=0.5,dehalo=0.2, S≈1.41523,D≈0.992366 and combined amplitude gain≈1.40442. This bin example treats weights as supplied; it does not claim every physical frequency can have both weights equal to 1.

## Placement

| Path | Composition |
|---|---|
| Effective one-frame uniform noise | Fuse Wiener×S×D using original R/q, then restore grid |
| Effective one-frame analytic/sampled noise | Wiener and grid restoration, then recompute enhancement g/q |
| T=2..5 | Temporal filtering, recover current spatial spectrum, then enhance |
| bt=−1 | Enhance original spatial spectrum without Wiener |
| bt=0 | Enhance an output copy of L; never feed it into recurrence |
| Effective preview | No enhancement gains |

If Wiener halves residual amplitude, power becomes one quarter. A subsequent S(q/4) generally differs from fused S(q). A profile path cannot be fused merely because its noise table happens to be constant. Preserve each branch's operation order; see [execution precision](../shared/execution-precision.md) and [Kalman](kalman.md).
