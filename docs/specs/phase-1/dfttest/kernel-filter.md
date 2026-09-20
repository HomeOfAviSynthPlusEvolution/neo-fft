# DFTTest spatial mean removal and frequency filters

Specification: DFT-FILTER-001. Evidence: static extraction from fixed sigma initialization, filter selection and scalar kernels.

## Inputs and outputs

Inputs: one block spectrum X; immutable G and wscale from [DFT-WIN-001](kernel-window-reconstruction.md); zmean; ftype in 0..4; finite nonnegative sigma,sigma2,pmin,pmax with pmin<=pmax; and finite f0beta>0. Output Y has identical logical compressed shape. Exact in-place mutation of X is allowed, while G and calibration are immutable and disjoint. Partial aliasing is forbidden.

All public numeric controls are converted to binary32 before calibration. No sigma curves, axis profiles or estimated noise are present in this phase.

## Per-bin constants

```text
wscalef = wscale if ftype<2, otherwise 1
A = sigma/wscalef
B = sigma2/wscalef
L = pmin/wscale
H = pmax/wscale
```

These are float divisions. A is a power threshold for types 0/1 and a gain parameter for types 2/3/4. The same public sigma is **not squared**. L,H are in the power units of the windowed transform. B is used only for type 3. Inactive constants are still validated for safe finite configuration.

## Mean removal and restoration

If zmean=true, compute g=real(X[0,0])/real(G[0,0]) before changing any bin, then `M[k]=g*G[k]` and `R[k]=X[k]-M[k]` for all bins. After the selected filter, restore `Y[k]=filtered(R[k])+M[k]`. With zmean=false, use R=X,M=0 and bypass the grid calculation.

This removes a DC-derived **window-shaped** spectrum, not only one coefficient and not the unweighted mean of input pixels. Different blocks have their own g. Retaining M for restoration must not accidentally overwrite the immutable G or another active block's removed component.

## Filter equations

For each bin let p=real(R)^2+imag(R)^2, calculated in float; epsilon is 1e-15f. Both components receive the same real gain.

| ftype | Rule |
| --- | --- |
| 0 | a=max((p-A)/(p+epsilon),0); gain selected by f0beta below; result=gain*R |
| 1 | If p<A, result=0; otherwise result=R. Equality retains the coefficient. |
| 2 | result=A*R |
| 3 | If L<=p<=H, result=A*R; otherwise result=B*R. Both boundaries are inclusive. |
| 4 | q=p+epsilon; gain=A*sqrt((q*H)/((q+L)*(q+H))); result=gain*R |

For type 0, evaluate float distance to the special exponents in this order:

```text
if abs(f0beta-1) < 0.00005f: gain=a
else if abs(f0beta-0.5f) < 0.00005f: gain=sqrt(a)
else: gain=pow(a,f0beta)
```

The comparisons are strict. At the exact boundary use the general branch unless float conversion changes the comparison. For other ftypes, f0beta has no numerical effect. It is still required to be finite and positive under this phase's safety contract.

Type 0 epsilon is only in the denominator; FFT3D places its epsilon in both numerator and denominator. Do not share a single Wiener kernel that erases this difference. Type 4 uses q including epsilon in its numerator and both denominator factors. Do not cancel/rearrange the rational expression or substitute an approximate reciprocal without validation. When H=0 and the expression remains representable the type-4 gain is zero; underflow/overflow yielding a non-finite expression is a controlled error.

## Examples

- R=3+4i gives p=25. With A=9, type 0 and exponent 1 give gain approximately 16/25=0.64, result 1.92+2.56i. Exponent 0.5 gives gain approximately 0.8 and result 2.4+3.2i.
- Type 1 with A=25 retains 3+4i; with A>25 it returns zero. This distinguishes `<` from `<=`.
- Type 3 with L=25,H=100,A=2,B=0.5 uses gain 2 at p=25 and p=100, and 0.5 outside the interval.
- Type 4 with p=25,L=0,H=100,A=1 gives gain approximately sqrt(0.8). Adding the FFT3D lower limit would be an error.
- For a constant windowed block X=(q/255)*G, zmean=true makes R=0 and restores X in exact arithmetic, even with type 2 sigma=0. With zmean=false,type 2 sigma=0 produces zero before synthesis.

## Errors and acceptance

Reject invalid domains, non-finite converted/calibrated constants and unusable G at creation. Content-dependent overflow/NaN is a frame error. Safe nonnegative domains and ordered pmin/pmax intentionally narrow loosely checked historical inputs; document negative/reversed cases as rejection changes, not equivalent outputs.

Test every ftype with zmean on/off, exponent 0.5/1/2 and neighboring values around both dispatch boundaries, exact and adjacent threshold powers, zero/all-DC input, signed/complex residuals, inactive padded bins and partial batches. Test threshold equality at kernel level with exact spectra in addition to end-to-end images; a VS round trip alone cannot reliably construct an exact floating spectral threshold.
