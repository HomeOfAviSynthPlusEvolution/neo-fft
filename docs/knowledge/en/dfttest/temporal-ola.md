# DFTTest temporal overlap-add

[Contents](../README.md)

`tmode=1` partitions time on a fixed block lattice. Each output receives contributions from all covering temporal blocks, independently of previous requests. Let T=tbsize,O=tosize,D=T−O and N be clip length.

[DFTTest overview](../dfttest.md) · [API](../../../api/en/dfttest.md)

## Lattice and endpoints

```text
A = -max(D,O)
s_j = A+j*D, j>=0
select every s_j<=n<s_j+T
target slice z = n-s_j
input slot k = clamp(s_j+k,0,N-1), k=0..T-1
```

Negative anchor A fixes the phase. Do not anchor it at the first requested frame or clamp s_j before forming slots. Source deduplication only reduces host fetches; repeated logical slots retain distinct window positions.

| T,O,n | Covering starts | Target z |
|---|---|---|
| 4,2,0 | −2,0 | 2,0 |
| 4,2,3 | 0,2 | 3,1 |
| 5,2,0 | −3,0 | 3,0 |
| 5,2,2 | 0 | 2 |
| 4,0,4 | 4 | 0 |

For N=5,T=4,O=2,n=4, starts 2 and 4 map to `[2,3,4,4]` and `[4,4,4,4]`. Repeated endpoints do not make these blocks identical. The ordinary union has at most min(N,2T−1) actual frames, plus active sampling dependencies; there is no Kalman-like persistent recurrence state.

## Temporal window and synthesis

Normalize raw window r by phase at step D: `tw[z]=r[z]/sqrt(sum r[k]², k≡z mod D)`. Sum z,z−D,... then z+D,... . O=0 still normalizes, turning nonzero r into its sign. Zero phase energy is an error.

Volume window `h=tw*sw_y*sw_x/sqrt(TB²)` determines the new wscale. After filtering the whole volume, its normalized inverse U contributes `(U[z]*TB²)*h[z]` at the target slice. Add in ascending start, spatial Y, spatial X order. smode=0 centers also add rather than replacing earlier temporal contributions. Convert/dither only once at the end.

## Two reproducible examples

For B=1,rectangular windows,T=4,O=2, tw=1/√2 throughout, h=1/√8,E=1/2,wscale=2. Two blocks each contribute half the original sample on an identity path; do not divide by 2 again.

For amplitudes `[0,2,4,6,8,10]`, zmean=true,ftype=2,sigma=0, output n=2 combines four-slot mean 3 from start 0 and mean 7 from start 2, each weighted by 1/2, yielding 5. The original sample was 4: the difference is mean filtering, not a synthesis gain error.

T=5,O=2 gives D=3 and tw=`[1/√2,1/√2,1,1/√2,1/√2]`, E=3/5,wscale=5/3. n=2 has one unit-weight contribution; n=3 has two half-weight contributions. A fixed overlap-count calibration would be wrong.

## Admission and precision

Require T in 1..15,T≤N,O in `[0,T)` and T divisible by D when overlap exceeds floor(T/2). T=6,O=4 is valid; T=5,O=3 is not. These checks also apply to `planes=[]`. tmode=0 still requires odd T, raw temporal windows and center output; the two modes are not guaranteed bitwise identical for every window.

One arithmetic configuration must support reverse, repeated and concurrent requests. No temporal-block cache currently exists; covering blocks are computed as defined. See [execution precision](../shared/execution-precision.md) for numerical and threshold boundaries.
