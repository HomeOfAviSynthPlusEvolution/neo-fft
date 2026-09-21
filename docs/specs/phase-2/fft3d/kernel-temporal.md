# FFT3D temporal DFT, degrid and PSD Wiener

Specification: F3D-TIME-001. Inputs are immutable phase-1 **raw windowed spatial spectra** X_j for T=btcur [source slots](../execution.md), phase-1 grid G and calibrated controls. Output is one spatial half-spectrum for n. Spatial shape is [bh,floor(bw/2)+1]; each bin's temporal transform is complex-to-complex, retaining **all T** temporal frequencies.

## Slot order and forward transform

| T | Chronological slots | Current index c |
| --- | --- | --- |
| 2 | prev1,cur | 1 |
| 3 | prev1,cur,next1 | 1 |
| 4 | prev2,prev1,cur,next1 | 2 |
| 5 | prev2,prev1,cur,next1,next2 | 2 |

At each spatial bin k, F[m,k]=sum(j=0..T-1) X_j[k]*exp(-2*pi*i*j*m/T). The forward transform is unnormalized; there is no temporal analysis window. The inverse must use c above; directly summing these chronological frequency bins reconstructs slot 0.

An equivalent unrolled implementation can use circular ordering A[r]=X[(c+r) mod T], so cur is at index 0. Its center inverse is a direct sum divided by T. Forward ordering and inverse phase must change together. The reference uses current-centered small transforms; account for frequency permutation/sign conventions in intermediate comparisons.

## Degrid uses current-frame DC, once per block

Use phase-1 spatial G and peak/sample units from [F3D-WIENER-001](../../phase-1/fft3d/kernel-wiener.md). If degrid!=0:

```text
g = (degrid*real(X_c[0,0]))/real(G[0,0])
M[k] = g*G[k]
gridT[k] = M[k]*float(T)
R[0,k] = F[0,k]-gridT[k]
R[m,k] = F[m,k]                         // m!=0
```

With degrid=0 bypass G/division and use M=0,R=F. The removed grid occupies **temporal frequency zero across every spatial bin**, not merely spatial DC. Its amplitude comes from the current block's DC, never neighborhood-average DC or independent g values for each source frame.

Subtracting the same M from all X_j is algebraically equivalent, but the compatibility baseline forms temporal sums first and subtracts T*M only from m=0; retain this order to limit drift. Independent per-frame degrid changes the filter when brightness varies. Cached X_j remain raw and immutable because one frame can be current in one request and a neighbor in another.

## Noise calibration and PSD

Let sigma_eff be the phase-1 sample-format-adjusted sigma (integer scaling 2^(bits-8), float scaling 1/255). With binary32 operations:

```text
norm = 1.0f/(float(bw)*float(bh))
noise = ((float(T)*sigma_eff)*sigma_eff)/norm
lower = (beta-1.0f)/beta
P = (real(R)*real(R) + imag(R)*imag(R)) + 1e-15f
gain = max((P-noise)/P,lower)
Rfiltered = gain*R
```

Apply to every temporal/spatial bin, with beta>=1. P is **squared magnitude plus epsilon**, never magnitude. Sigma is squared here, unlike DFTTest. Mathematically noise=T*bw*bh*sigma_eff^2, not spatial-window-energy calibration. The multiplication order follows make_shared_function_params; multiplying a pre-rounded 2D noise by T is not necessarily bit-identical. No half-spectrum doubling, extra T, upper gain clamp or division of PSD by volume is applied. Epsilon appears in numerator and denominator through P.

Validate finite noise/grid constants for every admitted effective T at creation, including fallback 1. Content-dependent non-finite PSD/transform/gain is a frame error. Process only logical bins, including tails; natural-alignment and padding rules remain in force.

## Inverse center and grid restoration

Restore Fout[0,k]=Rfiltered[0,k]+gridT[k], with other bins unchanged from Rfiltered, then:

```text
Y[k] = (1/T)*sum(m=0..T-1) Fout[m,k]*exp(+2*pi*i*c*m/T)
```

This restores M exactly once mathematically. In current-at-zero ordering, sum filtered bins plus gridT, then multiply by 1/T. Post-inverse +M is algebraically equivalent but may round differently. Follow with the phase-1 normalized **2D** inverse and spatial synthesis compensation bw*bh. The temporal inverse already supplied 1/T; synthesis must not multiply by T again. No temporal overlap-add/history exists.

T=1 boundary fallback runs the complete phase-1 path with T=1 noise, the same degrid, windows, conversion and properties. It is not endpoint replication, pass-through, a partially filled transform or a shorter T in 2..4.

## Diagnostic examples

- T=2, chronological x=[prev,cur]: F0=prev+cur,F1=prev-cur, target inverse=(Fout0-Fout1)/2. Rotated ordering gives F1=cur-prev and uses plus. An impulse only in prev must not become the target on a sigma=0 path.
- degrid=0,T=3,bw*bh=16,sigma_eff=2: noise=192. R=10+20i has P approximately 500 and beta=1 gain approximately 0.616. Magnitude-based filtering would incorrectly suppress it.
- T=3,G[k]=1,current DC ratio g=2, temporal values [1,2,6]: remove gridT=6 from F0=9, leaving R0=3. Nonzero temporal frequencies remain unchanged by grid subtraction; per-frame mean removal would erase temporal variation.
- Sigma=0 reconstructs cur within rounding for all T=2..5, degrid on/off. Use different complex values in every slot to expose reversed/rotated indexing.

See [acceptance](../acceptance.md). No patterns, sharpening or dehalo are opened by this contract.
