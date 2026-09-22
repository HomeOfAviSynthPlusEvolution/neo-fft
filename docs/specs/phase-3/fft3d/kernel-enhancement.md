# FFT3D: sharpen, dehalo and bt=-1

Specification: F3D-ENHANCE-003. Uses [interface](plugin.md) domains and phase-1 spatial reconstruction.

## Tables and gains

For each spatial half-spectrum bin use historical integer-half normalization, including odd-height behavior:

```text
dy = y if y<floor(H/2) else H-y
d2 = dy*dy*svr*svr/floor(H/2)^2 + x*x/floor(W/2)^2
Ws = 1-exp(-d2/(2*scutoff*scutoff))
rawWh = exp(-0.7*d2*hr*hr)-exp(-d2*hr*hr)
Wh = rawWh/max_yx(rawWh)
```

W,H>=2. At H=5, dy=[0,1,3,2,1]; do not substitute the estimator's symmetric coordinates. When dehalo=0 skip unused normalization; when active require finite positive maximum over logical bins. Both windows are zero at DC.

Let norm=1/(float(W)*float(H)), A=smin_eff^2/norm, B=smax_eff^2/norm, C=ht^2/norm. smin/smax use sigma's format factor. **ht does not**, matching pinned source; dehalo is not bit-depth-invariant. Constants and active tables must be finite. Calculations use binary32; evaluate exponentials into binary32 tables and preserve the written expression order except for independently validated SIMD rounding.

For spatial spectrum Z, compute g=degrid*Re(Z[0])/Re(G[0]) once per block (0 when disabled). R=Z-gG; q=Re(R)^2+Im(R)^2+1e-15:

```text
S(q) = 1 + sharpen*Ws*sqrt(q*B/((q+A)*(q+B)))
D(q) = (q+C)/((q+C)+dehalo*Wh*q)
Enhance(Z) = R*S(q)*D(q)+gG
```

Disabled gains are exactly 1 and their arithmetic is skipped. Both zero strengths bypass enhancement without subtract/add rounding. Require finite q and all used products, denominators, gains and results; do not suppress NaNs with min/max. B=0 gives S=1; ht=0 remains well-defined because q includes epsilon.

## Placement

| Path | Required composition |
| --- | --- |
| btcur=1, uniform | Fused: calculate R/q on original X, output R*Wiener(q)*S(q)*D(q)+gG |
| btcur=1, analytic/sampled pattern | Pattern Wiener/degrid, then Enhance on resulting spectrum, recomputing g/q |
| btcur=2..5 | Temporal Wiener/degrid, selective temporal inverse including 1/T, then Enhance on resulting spatial spectrum |
| bt=-1 | Enhance original spatial spectrum, no Wiener or temporal transform |
| effective pshow | Preview overrides these paths |

Do not enhance each temporal frequency independently or fuse the pattern path using original power. A second-stage degrid ratio comes from that stage's input DC. Model selection determines fused versus separate paths, not an inspection of whether P happens to be constant. All paths then use the existing normalized spatial inverse/synthesis.

bt=-1 with both strengths zero is the windowed reconstruction identity within existing rounding/quantization tolerances, independent of sigma. No input history is introduced.
