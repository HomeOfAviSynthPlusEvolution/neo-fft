# DFTTest UInt8 dither

Specification: DFT-DITHER-004. This operator consumes the fully reconstructed visible plane E in inherited 8-bit amplitude units. All E samples must be finite. dither=0 uses the phase-1 conversion; other sample formats also use their inherited conversion regardless of dither.

## Diffusion

For dither>=1 allocate two binary32 error rows current,next of visible width W. Initialize current to zeros before row 0; zero next before each row. Traverse y increasing, then x increasing on **every** row, not serpentine. Do not reset at a spatial block/strip boundary. For each pixel:

```
scale = float(dither - 1) + 0.5f
off = scale * 0.5f
v = (E[x,y] + current[x]) + 0.5f                     // dither == 1
v = (((E[x,y] + (u * scale)) - off) + current[x]) + 0.5f  // dither >= 2
D = clamp(trunc_toward_zero(v), 0, 255)
error = E[x,y] - float(D)
if x > 0:   next[x-1]    += error * 0.1875f
            next[x]     += error * 0.3125f
if x+1 < W: current[x+1] += error * 0.4375f
if x+1 < W: next[x+1]    += error * 0.0625f
```

For mode 1 evaluate `v=(E+current[x])+0.5f` directly, without adding an artificial zero noise term. For randomized modes do not regroup `E+u*scale-off` as `E+(u*scale-off)`: binary32 rounding differs. After completing a row swap current,next. Error outside the image is discarded; weights are not renormalized. Diffusion ends at the plane/frame boundary.

The historical error is **E-D**, not v-D and not E+incoming_error+noise-D. Use the clipped integer D in that subtraction. Preserve the stated update order and separate binary32 operations; do not contract/reassociate them. Clip before an unsafe float-to-integer cast: v<=0 yields 0, v>=255 yields 255, and only finite 0<v<255 is converted. This is equivalent to truncation followed by saturation for finite v, including negative values. Any non-finite active arithmetic fails the frame before publication. No integer conversion of NaN, infinity or an out-of-range float is permitted.

The final scan has a left-to-right and row-to-row dependency. No independent strip diffusion or per-worker error rows may replace it. Different selected planes can run independently; producing E in parallel remains allowed with the inherited deterministic accumulation order.

## Coordinate random sample

For dither>=2 use unsigned 32-bit arithmetic modulo 2^32, logical right shifts and these exact constants:

```
mix(v):
    v = v ^ (v >> 16)
    v = v * 0x7feb352d
    v = v ^ (v >> 15)
    v = v * 0x846ca68b
    return v ^ (v >> 16)

h = mix(uint32(seed) ^ 0xa511e9b3)
h = mix(h ^ uint32(n))
h = mix(h ^ uint32(p))
h = mix(h ^ uint32(y))
h = mix(h ^ uint32(x))
u = float(h >> 8) * 0x1p-24f
```

u is exactly representable binary32 in [0,1). Convert dither-1 only after the safe int32 subtraction (dither>=2); then apply the binary32 scale/off expressions above. No std::uniform_real_distribution, library-specific RNG, signed-overflow hash or hidden entropy is part of the contract. The hash is a deterministic noise source, not a cryptographic primitive.

Golden vectors, before multiplication by 2^-24:

| seed,n,p,y,x | h (hex) | h >> 8 |
| --- | --- | --- |
| 0,0,0,0,0 | 642a1d14 | 6564381 |
| 0,1,0,0,0 | d32e9dca | 13840029 |
| 1,0,0,0,0 | f383f4c5 | 15959028 |
| 17,23,2,5,7 | 9eb639f9 | 10401337 |
| 2147483647,2147483646,2,107,255 | 5064b117 | 5268657 |

Compare mode 1 directly with the pinned scalar reference on finite safe inputs. For modes >=2 compare against an independent implementation of this contract; the old workspace-seeded evolving RNG is intentionally not an exact-output oracle. Same-build output must be invariant under cache state, workers, selected-plane order, retries and random frame request order.
