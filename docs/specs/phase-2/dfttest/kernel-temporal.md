# DFTTest temporal window, mean, PSD and reconstruction

Specification: DFT-TIME-001. Inputs are an admitted [configuration](plugin.md), T ordered immutable [source slots](../execution.md), and phase-1 spatial geometry. Output is one spatial contribution per block for frame n. Let B=sbsize, T=tbsize, c=T/2 (integer), V=T*B*B, K=floor(B/2)+1. T is odd, including T=1 for regression.

## Window and calibration

Use the twelve raw window equations and I0 approximation in [DFT-WIN-001](../../phase-1/dfttest/kernel-window-reconstruction.md). They evaluate at j+0.5, not j or j/(L-1).

```text
tw[z] = raw_window(z,T,twin,tbeta)                // binary64; no temporal normalization
sw[j] = raw_window(j,B,swin,sbeta)                // binary64
if smode=1: normalize sw by phase-1 squared-overlap sums
h[z,y,x] = float(tw[z]*sw[y]*sw[x]*(1/sqrt(V)))    // double product; one storage cast
E = 0.0f
for z, then y, then x in ascending order:
    E = E + h[z,y,x]*h[z,y,x]                    // float multiply and accumulation
wscale = 1.0f/E
```

Require finite positive E and wscale. Energy uses the **stored float window**, excluding padding and inactive batches. Do not use energy of h*sqrt(V), invert twice, sum in double instead, or assume wscale=V. Analytically E=sum(tw^2)*sum(sw^2)^2/V; this is a check, not a replacement for the prescribed accumulation.

Sample conversion is unchanged: q=sample for UInt8; q=float(sample)*2^(-(bits-8)) for high-bit integer; q=sample*255 for Float32. Preserve signed float chroma; no midpoint subtraction. At each spatial origin gather **the same (x,y) block in all T slots**, then store a[z,y,x]=q*h[z,y,x]. Endpoint duplication repeats samples but each slot retains its own tw[z].

Forward [FFT-TIME-001](../kernel-fft.md) yields X of shape [T,B,K], scale 1. T is a transform axis, not an independent spatial-FFT batch axis.

## 3D zmean, before power evaluation

Prepare immutable G=FFT3D(255.0f*h) at creation using the same stored h, shape and backend as source blocks. For zmean=true require G finite and real(G[0,0,0]) finite and nonzero. A zero-DC window combination is a creation error; do not substitute epsilon or silently disable zmean. With zmean=false G and its division can be omitted.

For **each entire 3D block**, before overwriting its DC:

```text
g = real(X[0,0,0])/real(G[0,0,0])
M[k] = g*G[k]                        // every logical complex bin
R[k] = X[k]-M[k]
Y[k] = filter(R[k])+M[k]
```

With zmean=false use R=X,M=0. One g applies across all temporal/spatial frequencies in this block. Do not remove independent means per source frame, per temporal frequency, or only at one DC coefficient. This is a window-weighted DC ratio, not an unweighted pixel average. For constant q across the volume, X=(q/255)*G and R=0 in exact arithmetic. A time-varying sequence of spatial constants generally retains a temporal residual.

M survives filtering until restoration. Each active batch block needs its own g/M, or equivalent recomputation from retained g and immutable G. No partial aliasing of G, other active blocks or input/cache storage is permitted. Restoration precedes the 3D inverse and does not filter M.

## PSD and all five filters

Use binary32 public controls and the same [DFT-FILTER-001](../../phase-1/dfttest/kernel-filter.md) rules at every [kt,ky,kx] bin:

```text
wscalef = wscale if ftype<2 else 1.0f
A = sigma/wscalef
B2 = sigma2/wscalef
L = pmin/wscale
H = pmax/wscale
p = real(R)*real(R) + imag(R)*imag(R)
eps = 1e-15f
```

All constants/intermediates must be finite. PSD/p is **power**, not abs(R), sqrt(p), density per Hz or an average over T. Do not divide p by V or T, double compressed bins, or multiply thresholds by T again. DFTTest sigma is **not squared**; sigma2 independently defaults to 8.

| ftype | Filtered residual |
| --- | --- |
| 0 | a=max((p-A)/(p+eps),0); result a^f0beta * R with the exponent dispatch below |
| 1 | Zero if p<A, otherwise R; equality retains R |
| 2 | A*R; A is a gain, unscaled by window energy |
| 3 | A*R if L<=p<=H, otherwise B2*R; both bounds inclusive |
| 4 | qpower=p+eps; result A*sqrt((qpower*H)/((qpower+L)*(qpower+H)))*R |

Type 0 selects a when abs(f0beta-1)<0.00005f, sqrt(a) when abs(f0beta-0.5f)<0.00005f, otherwise pow(a,f0beta), in that order. Its epsilon is only in the denominator. Type 4 includes epsilon in numerator power and denominator factors. Preserve operation order and phase-1 overflow/error rules. No FFT3D beta floor is added.

## Normalized inverse and center synthesis

After restoring M, U=IFFT3D(Y) is normalized by **1/V exactly once**. Extract temporal slice c, which always maps to frame n, including boundaries. Use float operations in this order:

```text
contribution[y,x] = (U[c,y,x]*float(V))*h[c,y,x]
```

V compensates the legacy unnormalized inverse; omitting V loses volume gain and using only B*B loses a factor T. There is no extra 1/T or inverse-window division. A selective-center inverse is allowed only if equivalent to this full definition within the frozen budget.

For smode=1 add all B*B contributions into the zeroed **2D** spatial accumulator in phase-1 origin Y-then-X order, then crop. For smode=0 emit only local (B/2,B/2) at its corresponding output location. Retain phase-1 padding, output conversion and safe quantization. No contribution is published to another output frame.

On an identity spectrum, ideal overlap gain is tw[c]^2; for smode=0 it is tw[c]^2*sw[B/2]^4. With rounded tables the authoritative overlap gain is the sum of V*h[c,local_y,local_x]^2 over covering blocks. Do not force unity. Odd T puts r at 1/2, but **window 6 has tw[c]=0.2810639+0.5208972+0.1980399=1.000001**, giving ideal temporal gain 1.000002000001. Window 9 remains asymmetric away from center. Window 5 also uses its actual coefficient sum, not a hardcoded center 1.

## Hand-checkable cases

- B=2,O=1,smode=1,swin=twin=7,T=3: sw=1/sqrt(2), h=1/(4*sqrt(3)), ideal E=1/4,wscale=4. sigma=8,ftype=0 gives A approximately 2. Synthesis uses V=12 and four contributions q/4.
- Same spatial setup,T=3,twin=0: tw=[1/4,1,1/4], ideal E=3/32,wscale=32/3,A=3/4 for sigma=8. Center identity gain remains 1; this distinguishes full-volume from center-only noise calibration.
- B=1,smode=0,swin=twin=7,T=3,q=[1,2,6]: h=1/sqrt(3). With zmean=true,ftype=2,sigma=0 the target reconstructs 3 (temporal mean), not original value 2. With zmean=false it reconstructs 0. Float fixtures supply q/255.
- T=1 reduces to phase 1, including all twelve twin choices and actual center gains; it must not bypass window evaluation.

See [acceptance](../acceptance.md). Analytic examples allow float rounding and do not promise bit identity between FFT backends.
