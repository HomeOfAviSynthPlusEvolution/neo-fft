# DFTTest spectral gains, curves and sampling

[Contents](../README.md)

Let X be the windowed space/time spectrum and G=FFT(255h). With zmean, `g=Re(X[0])/Re(G[0])`, `M=gG`, `R=X−M`; otherwise M=0,R=X. Filtering acts on R and restores M afterward.

[DFTTest overview](../dfttest.md) · [Complete parameters](../../../api/en/dfttest.md)

## Five spectral rules

Set `p=Re(R)²+Im(R)²`, ε=1e−15, wscale=1/Σh². For primary model v, A=v/wscale for types 0/1 and A=v for 2/3/4. B applies the same divisor to sigma2. L=pmin/wscale,H=pmax/wscale.

| ftype | Residual gain |
|---|---|
| 0 | a=max((p−A)/(p+ε),0), gain=a^f0beta |
| 1 | 0 when p<A, otherwise 1; equality retains |
| 2 | A |
| 3 | A when L≤p≤H, otherwise B; both endpoints included |
| 4 | q=p+ε, gain=A√(qH/((q+L)(q+H))) |

Type 0 uses a directly when float `abs(f0beta−1)<0.00005`; otherwise √a when `abs(f0beta−0.5)<0.00005`; otherwise pow. Comparisons are strict. Epsilon placement differs from FFT3D, so a shared algebraically simplified Wiener expression would be wrong.

R=3+4i,p=25,A=9 gives type-0 gain≈0.64 and result≈1.92+2.56i for exponent 1, or gain≈0.8 for exponent 0.5. Type 1 retains R at A=25 and clears only for A>25. Type 3 with L=25,H=100,A=2,B=0.5 uses 2 at p=25.

## Curves to v

All-empty curves use sigma directly. Otherwise d=(T>1?1:0)+(B>1?2:0). For axis length L>1, `f(i,L)=min(i,L−i)/floor(L/2)`; length 1 uses 0. T=3 maps to `[0,1,1]`.

Raise curve values to exponent e **before** linear interpolation. Shared slocation uses e=1/d for ssystem=0, e=1 for ssystem=1. Axis-only curves use e=1/d for both systems. Empty axes use sigma^e; exact knots return their stored values.

ssystem=0 multiplies temporal/Y/X interpolants; singleton axes contribute 1. ssystem=1 computes `rho=sqrt((ft²+fy²+fx²)/d)` and reads the temporal table at rho. With axis-only curves, ssx/ssy do not supply radial ordinates, so use slocation for ordinary radial profiles. At d=0, take slocation DC, else sst DC, else sigma, with exponent 1.

For T=1,B=3,d=2, shared separable curve `[0,4,1,16]` first becomes `[0,2,1,4]`. At fx=1,fy=0, v=4×2=8. Do not interpolate untransformed values and then apply an arbitrary root. Shared radial mode instead interpolates 4 to 16 at rho=1/√2.

## Sampling replaces the primary model

For types 0/1, nonempty nlocation overrides primary sigma/curves. Sampling h2 uses raw spatial and temporal windows without either overlap normalization, retaining 1/√(TB²). Each tuple consumes an original in-bounds rectangle over T consecutive frames in DFTTest amplitude units. With zmean, remove its window-shaped mean using G2=FFT(255h2).

```text
wscale2 = 1/sum_float(h2*h2)
Psum[k] = sum_in_tuple_order(abs(R_sample[k])**2)
A[k] = Psum[k] * ((1/M)*(wscale2/wscale)*alpha)
```

M counts repeated tuples. A is already calibrated power: no second division by wscale, squaring or T multiplier. It is shared across output planes and does not replace sigma2/pmin/pmax. Two powers 4 and 8, M=2, energy ratio 0.5 and alpha=5 give A=15.

Overridden raw arrays still validate, but overridden derived curves need not be built. Negative values, duplicate float32 knots, missing endpoints, out-of-bounds patches or non-finite active powers fail. Sampling and filtering means use their own templates. Results follow [execution precision](../shared/execution-precision.md).
