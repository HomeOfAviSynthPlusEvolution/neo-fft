# UInt8 quantization and deterministic dither

[Contents](../README.md)

DFTTest enables dither only on selected UInt8 planes, after complete space/time synthesis and crop. It does not change spectra or windows. Higher integer depths and float keep ordinary conversion for every valid mode; negative modes/seeds fail even when inactive.

[API](../../../api/en/dfttest.md) · [Sample units](../shared/sample-domain.md)

## Scan and error

Clear two error rows for each plane/frame. Scan every row left to right, rows top to bottom: no serpentine order or resets at block boundaries. Let E be reconstructed 8-bit amplitude and e incoming error:

```text
dither=0: v = E+0.5
dither=1: v = (E+e)+0.5
dither>=2:
    scale = float(dither-1)+0.5
    off = scale*0.5
    v = (((E+u*scale)-off)+e)+0.5
D = clip(trunc(v),0,255)
error = E-D
```

Error is original E minus clipped D, **not** noisy/error-adjusted v−D. Update next-row left, next-row same column, current-row right and next-row right with weights 3/16,5/16,7/16,1/16 respectively. Discard out-of-image contributions without renormalization. Frames are independent.

For two pixels E=[10.4,10.4] and dither=1, the first yields D=10,error=0.4 and sends 0.175 rightward. The second has v=10.4+0.175+0.5=11.075,D=11. Its own error is −0.6, not 0.075.

## Coordinate-derived noise

n is output frame, p the actual plane index and x/y visible plane coordinates. They do not depend on planes-array order, block number or worker identity.

```text
mix(v):                              # uint32, modulo 2^32
    v = (v ^ (v >> 16)) * 0x7feb352d
    v = (v ^ (v >> 15)) * 0x846ca68b
    return v ^ (v >> 16)
h = mix(uint32(seed) ^ 0xa511e9b3)
for coordinate in [n,p,y,x]: h = mix(h ^ uint32(coordinate))
u = float(h >> 8) * 2^-24
```


At seed=n=p=y=x=0, h=0x642a1d14,h>>8=6564381,u=6564381/16777216. For dither=2,scale=1.5,off=0.75, so the perturbation is `1.5u−0.75`, approximately −0.1631. Requesting another frame first cannot change this u.

Modes 0/1 do not use seed or advance hidden random state. Modes ≥2 use this project's deterministic contract, not an exact reproduction of the old workspace RNG. Floating operation order can change half-LSB decisions: do not replace `E+u*scale-off` with `E+(u*scale-off)`. Noise can be generated with SIMD, but diffusion is dependent and currently scans sequentially on the caller.

Finiteness checks and safe saturation precede unsafe integer casts. Failure does not publish partial frames. Repeated, reverse and concurrent requests in one arithmetic configuration must agree; see [execution precision](../shared/execution-precision.md).
