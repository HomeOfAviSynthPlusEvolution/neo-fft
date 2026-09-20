# DFTTest windows, calibration and reconstruction

Specification: DFT-WIN-001. Evidence: static extraction from window construction, engine calibration and scalar spatial processing. The normalized-inverse compensation is a derived requirement.

## Inputs and outputs

Inputs: admitted B,O,smode; swin,twin in 0..11; nonnegative finite sbeta,tbeta; T=1; sample format. Outputs: a binary32 2D window h, its energy calibration, analysis blocks and quantized/cropped output. Input/table/output ownership follows [MEM-001](../kernel-plane-geometry.md); accumulation scratch is exclusive and initially zero.

## Raw window

At index j of length L, let n=j+0.5, r=n/L and theta=2*pi*r. Calculate in binary64. Cosine terms below are cos(m*theta).

| ID | Name in reference | Raw w(j,L,beta) |
| --- | --- | --- |
| 0 | Hann | 0.5 - 0.5*cos(theta) |
| 1 | Hamming | 0.53836 - 0.46164*cos(theta) |
| 2 | Blackman | 0.42 - 0.5*cos(theta) + 0.08*cos(2*theta) |
| 3 | Blackman-Harris, 4 term | 0.35875 - 0.48829*cos(theta) + 0.14128*cos(2*theta) - 0.01168*cos(3*theta) |
| 4 | Kaiser-Bessel | I0a(pi*beta*sqrt(1-(2*r-1)^2)) / I0a(pi*beta) |
| 5 | Blackman-Harris, 7 term | a0 + a1*cos(theta) + ... + a6*cos(6*theta), coefficients below |
| 6 | Flat top | 0.2810639 - 0.5208972*cos(theta) + 0.1980399*cos(2*theta) |
| 7 | Rectangular | 1 |
| 8 | Bartlett | (2/L)*(L/2-abs(n-L/2)) |
| 9 | Bartlett-Hann label | 0.62 - 0.48*(r-0.5) - 0.38*cos(theta) |
| 10 | Nuttall | 0.355768 - 0.487396*cos(theta) + 0.144232*cos(2*theta) - 0.012604*cos(3*theta) |
| 11 | Blackman-Nuttall | 0.3635819 - 0.4891775*cos(theta) + 0.1365995*cos(2*theta) - 0.0106411*cos(3*theta) |

For ID 5, `(a0,...,a6)` is:

```text
(0.27105140069342415,
 -0.433297939234486060,
 0.218122999543110620,
 -0.065925446388030898,
 0.010811742098372268,
 -7.7658482522509342e-4,
 1.3887217350903198e-5)
```

I0a is the reference's finite approximation, not an arbitrary library Bessel function. Start numerator=1, denominator=1, total=1. For k=1..14, multiply numerator by argument/2 and denominator by k, set v=numerator/denominator and add v*v to total. Stop after that addition when v<=1e-8; otherwise continue through k=14. Return total. Nonnegative beta avoids the reference's signed early-stop anomaly. Reject overflow/non-finite window calculations at creation.

ID 9 deliberately has **no absolute value** around r-0.5. It is asymmetric; do not replace it with the conventional named window or force spatial symmetry. Window names are descriptive labels, while these equations define compatibility. Negative finite raw coefficients are allowed.

## Normalization and energy

Set raw[j]=w(j,B,sbeta). For smode=0, sw=raw. For smode=1 with S=B-O:

```text
D[j] = sum(raw[k]^2 for 0<=k<B and k congruent to j modulo S)
sw[j] = raw[j]/sqrt(D[j])
```

For a reproducible scalar table, sum in the reference order j,j-S,... down to 0, followed by j+S,j+2S,... below B. Require every D[j] finite and positive. There is no spatial normalization in center-sample mode.

Set tw=w(0,1,tbeta) using twin; T=1 does not automatically override twin to rectangular. There is no temporal overlap normalization with tmode=0. Set N=B^2 and store:

```text
h[y,x] = float(tw*sw[y]*sw[x]/sqrt(N))
E = sum_float_row_major(h[y,x]*h[y,x])
wscale = float(1/E)
```

Require finite positive E and wscale. Store h only after computing its products in double. E uses stored float h and float sequential accumulation, as in the reference. Power thresholds use division by wscale, not an assumed ideal analytical E; see [filter](kernel-filter.md).

## Analysis units

Integer 8-bit input q=sample. High-bit integer q=float(sample)*2^(-(bits-8)). Float32 input q=sample*255. Analysis is q*h with the conversion multiplication performed first. No YUV chroma midpoint is subtracted, and no legal-range remapping is performed. Signed float chroma therefore enters signed and scaled by 255.

Prepare the immutable mean template `G=FFT(255*h)` using stored h. G uses the same FFT backend/shape as source blocks. If zmean is enabled, require finite nonzero real G[0,0]. This template is used for every plane, regardless of bit depth or YUV/RGB identity.

Window/mode combinations that have zero template DC are outside the zmean domain even if the window ID itself is supported. For example, flat-top B=4 with O=0 in smode=1 normalizes to spatial signs [-1,1,1,-1], so the template has zero DC. Reject zmean=true at creation rather than divide by zero; zmean=false can still use that finite window.

## Synthesis and output

Let U be the **normalized** inverse transform. Each local contribution is `(U*N)*h`, with float operations in that order. N compensates the reference's unnormalized inverse. It is not a second inverse normalization. Keeping the analysis h unchanged preserves the spectral epsilon and power calibration.

In smode=1 add contributions for all local samples in origin Y-then-X order, then crop. The exact identity gain on the fully covered crop is tw^2 because each spatial phase has sum(sw^2)=1. Do not divide by this gain. In smode=0 emit only local (c,c), c=floor(B/2); its exact identity gain is `N*h[c,c]^2`, or `tw^2*sw[c]^4` before table rounding. No corrective center normalization is applied.

Let z be the final cropped accumulator in 8-bit amplitude units. With dither=0:

| Output format | Conversion |
| --- | --- |
| UInt8 | clip(trunc(z+0.5),0,255) |
| UInt10/12/14/16 | clip(trunc(z*2^(bits-8)+0.5),0,2^bits-1) |
| Float32 | z*float(1/255), without clipping |

Use safe saturation before any out-of-range integer cast. Float output preserves negative/out-of-nominal values. A copied plane undergoes none of these conversions. Selected non-finite samples or non-finite content-dependent calculations fail the frame.

## Examples and acceptance

- B=2,O=1,smode=1,swin=twin=7: sw=1/sqrt(2),h=1/4,E=1/4,wscale=4. Four identity contributions q/4 reconstruct q. Public sigma=8 means a type-0 threshold of 2, not 8^2 and not 8.
- At B=2, ID 9 produces raw [0.74,0.50]; these values distinguish the source equation from the conventional absolute-value window. It must not be mirrored into a symmetric table.
- A 10-bit sample 512 enters as 128. Float chroma -0.25 enters as -63.75 and returns -0.25 on an identity path, subject to float rounding, without clipping to zero.
- B=1,smode=0,swin=7,twin=7 gives h=1,E=1 and a direct per-pixel spectral operator. ftype=2,sigma=1,zmean=false is identity.

Verify all twelve spatial and T=1 temporal window choices, Kaiser beta including zero, asymmetric window 9, zero overlap and legal heavy overlap, float-table energy, midpoint scale and output half-LSB boundaries. Identity tests must use the actual window gain; a non-unit center-window gain is not automatically a filtering bug.
