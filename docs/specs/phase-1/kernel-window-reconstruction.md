# Window application and reconstruction primitives

Specification: WIN-001. Evidence: shared operator design; algorithm-specific normalization remains in the two domain specs.

## Inputs and outputs

Inputs are checked source/block planes, finite binary32 coefficient tables and explicit block origins. Analysis produces a disjoint float block; synthesis either adds a full block to a float accumulator or emits its specified center. Tables and source planes are immutable. The accumulator is initialized to zero before its first contribution and belongs to one request.

For source conversion q, analysis a, synthesis s and normalized inverse U:

```text
B[y,x] = q(source[origin_y+y,origin_x+x]) * a[y,x]
U = inverse(filter(forward(B)))
accumulator[origin_y+y,origin_x+x] += U[y,x] * s[y,x]
```

Center sampling instead emits exactly `U[cy,cx]*s[cy,cx]`; it does not overlap-add all pixels and then divide by the number of contributing blocks.

## Algorithm calibration

| Rule | FFT3D | DFTTest, T=1 |
| --- | --- | --- |
| Working sample | Native sample minus integer YUV chroma midpoint | 8-bit amplitude scale, no centering |
| Analysis | A=Ax*Ay | h=tw*sw_y*sw_x/sqrt(N) |
| Synthesis with normalized IFFT | Sx*Sy | N*h |
| Noise calibration | sigma squared times N | ftype-dependent sigma and window-energy scaling |
| Output rule | Integer round/clip; float clip to [0,1] | Integer round/clip; float scale back without clip |

Do not unify these units because both pipelines use float storage. In particular DFTTest's h includes 1/sqrt(N); moving that factor out of analysis changes the meaning of its fixed spectral epsilon.

## Coverage and operation order

With filtering replaced by identity, full-block reconstruction has weight `w[p]=sum_blocks(a[p]*s[p])`. This is a proof/test quantity, **not** a request to divide every output by w. FFT3D's visible crop has w=1. DFTTest overlap-add has w=tw^2 in exact arithmetic; for default rectangular temporal window at T=1 it is 1. Center sampling has its own window gain.

FFT3D combines horizontal contributions before vertical contributions. DFTTest visits origins in increasing Y then X and accumulates their contributions in that order. Batch boundaries do not reorder contributions. Optimized paths may only change rounding within the accepted numerical budget; an arbitrary parallel reduction is not part of the scalar baseline.

Store coefficients and accumulators as binary32. Algorithm specs override the general binary64 coefficient-generation default where the fixed source uses float calculations. Quantize once after all contributions; never quantize a block before adding it. No clipping occurs between FFT, filtering and reconstruction unless explicitly specified.

## Examples and failures

Two overlapping contributions with analysis/synthesis pairs (a,s)=(1/2,1) and (1/2,1) reproduce a constant q since q/2+q/2=q. Two pairs (1/2,1/2) yield q/2, and require a different window design; silently dividing by observed weight would conceal that error.

For DFTTest B=2, rectangular spatial and temporal windows, O=1, normalized sw=1/sqrt(2) and h=1/4. N*h=1. Each identity block contributes q/4; four blocks reproduce q. Omitting N from synthesis produces q/4 instead.

Reject non-finite/zero window denominators or invalid coverage at plan creation. Content-dependent non-finite accumulation fails the frame. All visible output pixels must be initialized, including unprocessed planes. Window and grid tables cannot alias a mutable accumulator.
