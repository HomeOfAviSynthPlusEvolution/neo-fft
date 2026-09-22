# DFTTest: curve and sampled sigma tables

Specification: DFT-PROFILE-003. Shape is [T,S,floor(S/2)+1], x contiguous, one real primary value per complex bin, shared by real/imaginary components. Existing h, wscale=1/sum_float(h*h), G=FFT(255*h), zmean and ftype formulas remain as in [phase 2](../../phase-2/dfttest/kernel-temporal.md).

## Curves

Let d=(T>1 ? 1:0)+(S>1 ? 2:0); a=1/d when d>0. Raise each curve ordinate v to its chosen exponent e **before interpolation**. Empty curves supply endpoints (0,sigma^e),(1,sigma^e). Between sorted knots (u,v),(w,r), evaluate v*(1-t)+r*t with t=(f-u)/(w-u). Exact knots return their stored ordinate. Do not interpolate raw values then take a root.

For length L>1 use f(i,L)=min(i,L-i)/float(floor(L/2)); for L=1 set f=0 and the separable axis factor to 1. T=3 therefore gives [0,1,1], not [0,2/3,2/3].

If all curve arrays are empty, use scalar sigma without entering the curve path, irrespective of ssystem. Otherwise:

| Selection | Prepared tables |
| --- | --- |
| slocation nonempty | All axes use slocation, exponent a for ssystem=0 or 1 for ssystem=1 |
| slocation empty | T/Y/X use sst/ssy/ssx, exponent a for **both** systems; empty axes use default sigma |

For d>0:

- ssystem=0: V=It(ft)*Iy(fy)*Ix(fx), singleton axes contribute 1.
- ssystem=1: rho=sqrt((ft*ft+fy*fy+fx*fx)/d), V=It(rho). Axis-only mode uses the temporal table even when T=1; ssx/ssy do not supply radial ordinates. This unusual pinned behavior must be tested. For a conventional radial curve, supply slocation.

For d=0 (T=S=1) use exponent 1 and one DC ordinate: slocation at 0 if present; otherwise sst at 0 if present; otherwise sigma. ssx/ssy alone do not alter this DC value. This is an explicit safe singleton extension, not reference equivalence.

Primary A[k]=V[k]/wscale for ftype=0/1 and V[k] for ftype=2/3/4. sigma2 uses the same type-dependent scalar divisor; pmin/pmax always divide by wscale. Apply existing ftype equations per bin. Do not square V or add a volume/T multiplier. For the consumed profile, check finite pow/interpolation/products/divisions and nonnegative tables at construction, treating negative zero as zero. A profile overridden by active nlocation is not constructed; its raw arrays are still validated under [execution](../execution.md). Use binary32 table/interpolation arithmetic and the existing tolerance policy for libm/SIMD rounding.

## Noise samples and energy calibration

Only active for ftype=0/1 with nonempty nlocation; M counts tuples including repetitions. Build a second analysis window h2 with configured swin/twin/sbeta/tbeta and half-sample evaluations, with **both spatial and temporal overlap normalization disabled**, regardless of output smode/sosize. It still includes 1/sqrt(T*S*S). Do not change output filtering window h.

```text
E2 = row-major binary32 sum of h2*h2
wscale2 = 1/E2
G2 = FFT(255*h2)                        # full T*S*S transform
```

Require positive finite E2 and finite h2/G2; zmean=true also requires nonzero Re(G2[0]). For tuple (fn,p,y,x), gather its in-bounds S*S plane-p patch from frames fn..fn+T-1, in chronological order. Convert to normal DFTTest 8-bit amplitude units (integer divide by 2^(bits-8), float multiply by 255), with no FFT3D midpoint subtraction, then multiply by h2. Let Xj be the forward FFT. If zmean=true use Rj=Xj-Re(Xj[0])/Re(G2[0])*G2; otherwise Rj=Xj. Capture DC before any in-place modification.

```text
Psum[k] = sum_j (Re(Rj[k])^2+Im(Rj[k])^2)
A[k] = Psum[k] * ((1/float(M)) * (wscale2/wscale) * alpha)
```

Accumulate tuples in supplied order, bins in canonical order, in binary32. The energy ratio converts sampling-window FFT power to filtering-window power. A is already calibrated for ftype=0/1: **no second division by wscale**, sqrt, squaring, T factor or half-spectrum multiplicity. This table replaces primary A on all selected output planes. Normal filtering continues to use G/h and its normal mean restoration; G2/h2 only estimate noise.

Check finite sampled values and every accumulated/derived intermediate. zmean=false intentionally includes DC power. Constant samples with zmean=true yield zero residual in exact arithmetic, subject to established FFT rounding. Duplicating the entire tuple list preserves the average within rounding tolerance.
