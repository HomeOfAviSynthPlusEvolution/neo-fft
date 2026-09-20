# FFT3D windows, sample conversion and reconstruction

Specification: F3D-WIN-001. Evidence: static extraction from overlap-window generation and cover/overlap transforms; normalized-IFFT compensation is derived.

## Inputs and outputs

Inputs: the admitted [geometry](kernel-geometry.md), immutable source samples, format/plane identity, wintype in {0,1,2}, and transformed blocks. Output: window tables, analysis blocks and a fully initialized visible destination plane. Source, destination and mutable scratch do not alias. Tables are read-only; spectral filtering may operate in place as specified separately.

## Axis windows

For axis length B and overlap O, use 1 in the non-overlap interior. If O=0, the entire analysis/synthesis window is 1. Otherwise, define overlap values for j=0..O-1:

```text
cL[j] = cos(pi*(j-O+0.5)/(2*O))
cR[j] = cos(pi*(j+0.5)/(2*O))
```

At positions j in the left overlap and B-O+j in the right overlap choose c=cL[j] and cR[j] respectively:

| wintype | Analysis a | Synthesis s |
| --- | --- | --- |
| 0 | c | c |
| 1 | sqrt(c) | a*a*a |
| 2 | 1 | c*c |

Use binary32 pi, intermediates, cosine/square-root results and table storage, following the fixed source. Calculate type-1 synthesis from its stored analysis value; do not replace it with an independently evaluated power. A precision-changing implementation requires differential evidence, not an assumption that more precision is interchangeable.

For 2D local coordinates, A[y,x]=ay[y]*ax[x], S[y,x]=sy[y]*sx[x]. Every overlapping axis pair has `aL*sL+aR*sR=1` in exact arithmetic because cL^2+cR^2=1. With O<=B/2 there are at most four contributors to any visible pixel.

## Sample and noise units

Integer samples use their native numeric value. Subtract `base=2^(bits-1)` only for integer YUV U/V planes; all other planes, including float YUV chroma and RGB, use base=0. The new scalar analysis stores `((float(sample)-base)*ay)*ax` using the plane's proper windows. This fixes one order across all block positions; the reference's corner branches sometimes multiply the windows first, so rounding differences must satisfy the numerical limits in the [acceptance protocol](../acceptance.md#numerical-acceptance). There is no limited-range expansion, division by peak or 255 multiplier for float input.

The public sigma unit is 8-bit standard deviation. Effective sigma is `float(sigma)*factor`, with factor 1 for 8-bit, 2^(bits-8) for high-bit integer, and float(1/255) for float32. This scaling is per plan/plane; never multiply a shared configuration repeatedly while creating three plane engines. Spectral power calibration is defined in [Wiener](kernel-wiener.md).

## Reconstruction

Let U be the normalized inverse transform of each filtered block. For each visible cover position, sum `U[y,x]*S[y,x]` over covering blocks. Combine left/right contributions first, then top/bottom contributions. No extra 1/N and no run-time division by a weight image are applied. The old exterior-cover shortcuts do not affect the admitted visible crop.

Use float arithmetic throughout reconstruction. With identity filtering the crop reconstructs sample-base up to numerical error. Add the midpoint back only at final conversion.

For integer output with peak=2^bits-1, form `v=(reconstructed+0.5)+base`, truncate toward zero, and clamp to [0,peak]. Implement saturation before an unsafe float-to-integer conversion when v exceeds the representable range; for finite v this must produce the same clipped result. The order of +0.5 and +base is specified because the reference performs it in float.

For float output, return `clamp(reconstructed,0,1)` with no midpoint or rounder. **This includes float YUV chroma.** Negative chroma and out-of-range float values are clipped by the fixed reference; this phase preserves that behavior. DFTTest has a different float output rule. Unselected planes are copied bitwise and are never clipped.

Reject non-finite selected input samples or intermediate/output calculations with a frame error, not a NaN-dependent integer cast. Values outside nominal float range are otherwise valid inputs and undergo the stated final clipping.

## Examples

- B=4,O=2: cL=[sin(pi/8),cos(pi/8)] and cR=[cos(pi/8),sin(pi/8)]. At the first overlap position the reconstruction weights sum to approximately 0.1464466+0.8535534=1 for all three window types.
- A neutral 10-bit chroma value 512 enters as 0; sigma=2 becomes 8. Identity reconstruction returns 512. The same numeric sample in an RGB plane enters as 512.
- A float chroma input -0.25, preserved by identity filtering, becomes 0 on a selected plane; a copied plane retains -0.25. A value 1.25 becomes 1.
- Integer luma reconstructed as 12.49 yields 12, and 12.5 yields 13. There is no ties-to-even pixel quantizer.

## Acceptance

Test each window type with O=0, asymmetric overlaps and half overlap; verify the visible coverage weight, impulse reconstruction and constant planes independently of spectral denoising. Include integer chroma midpoints, RGB, float signed chroma, clipping limits and half-LSB boundaries. Compare actual table precision and final output separately so a window error cannot hide behind quantization.
