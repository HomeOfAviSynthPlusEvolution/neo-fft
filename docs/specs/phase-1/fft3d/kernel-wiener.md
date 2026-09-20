# FFT3D spatial Wiener and degrid

Specification: F3D-WIENER-001. Evidence: static extraction from the fixed 2D scalar and Highway paths. Temporal, pattern, sharpening and dehalo operations are excluded.

## Inputs and outputs

Inputs: one windowed block spectrum X, Bx,By, effective sigma, beta>=1, degrid>=0, and the immutable grid spectrum G when degrid is nonzero. Values must be finite and representable in binary32. Output Y has the same compressed shape. Exact X-to-Y in-place filtering is allowed; partial overlap or aliasing G is not. Process only logical bins.

Let N=Bx*By and compute `norm=1.0f/(float(Bx)*float(By))`. Preserve reference calibration order: `noise=(sigma_eff*sigma_eff)/norm`. Its mathematical value is N*sigma_eff^2. Do not replace N with window energy, apply Hermitian doubling, or square the public sigma twice.

## Grid component

Generate G by forwarding a block of the analysis window multiplied by peak, with peak 2^bits-1 for integer formats and 1 for float. This template is generated **without chroma midpoint subtraction**, even when processing chroma. Use the same analysis multiplication order and FFT backend as source blocks. Its real DC must be finite and nonzero.

For each block, before changing X[0,0]:

```text
g = degrid * real(X[0,0]) / real(G[0,0])
M[k] = g * G[k]
R[k] = X[k] - M[k]
```

For degrid=0 bypass the grid calculation, with R=X and M=0. Default degrid=1 removes/restores the full DC-derived window-shaped component. This is not subtraction of only the DC coefficient and not an arithmetic spatial mean of unwindowed pixels.

## Wiener rule

For every stored bin, using binary32 arithmetic:

```text
P = (real(R)*real(R) + imag(R)*imag(R)) + 1e-15f
lower = (beta-1)/beta
gain = max((P-noise)/P, lower)
Y = gain*R + M
```

No separate upper clamp is needed in the admitted nonnegative-noise domain. The epsilon occurs in P in both numerator and denominator. With beta=1, gain is zero when P<=noise; beta>1 retains a floor. Restoring M happens after filtering and before inverse FFT. Keep scalar divisions exact to the chosen float operation; reciprocal approximations/FMA require their own verification.

The gain is real and preserves the input's Hermitian compatibility. Never multiply spectral DC/Nyquist bins by an additional normalization constant. Degrid coefficients are computed independently per block; a preceding block's DC cannot be reused.

## Examples

- N=16, 8-bit sigma=2 gives noise=64. R=10+0i gives P approximately 100 and beta=1 gain approximately 0.36, producing 3.6. beta=2 gives gain 0.5 and output 5 instead.
- A bin R=3+4i has power approximately 25. With noise=25,beta=1 it is suppressed (apart from representable epsilon effects). With sigma=0 the computed ratio P/P is 1 even at a zero bin.
- If a constant plane makes X=(c/peak)*G and degrid=1, R=0 in exact arithmetic and Y=X for any admitted sigma. degrid=0 generally attenuates it. This checks subtraction/restoration across all bins, not just DC.
- With 10-bit sigma=2, effective sigma=8 and N=16, noise=1024. For float sigma=2, effective sigma is about 0.00784314 and noise about 0.000984237.

## Errors and acceptance

Negative/non-finite sigma or degrid, beta<1, unusable grid DC and non-finite derived noise are creation errors. Content-dependent overflow/non-finite powers are frame errors. These safety requirements are deliberate validation tightening; the old parser's acceptance of arbitrary numeric values does not define safe output for them.

Test gain-floor and threshold branches, zero residual, degrid 0/0.5/1, different block DC values in one batch, odd/even and non-square spectra, and scalar/Highway tails. Test the complete path with each pixel format; a correct isolated Wiener formula does not establish correct sigma scaling or overlap reconstruction.
