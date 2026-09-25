# Windows, energy and reconstruction

[Contents](../README.md)

Analysis windows control block edges entering the FFT; synthesis windows weight inverse contributions. Treat them as a pair. Dividing by the number of contributors afterward can conceal a scale error.

## FFT3D shoulder windows

For overlap O and j=0..O−1, `cL=cos(π(j−O+0.5)/(2O))`, `cR=cos(π(j+0.5)/(2O))`. The interior is 1; O=0 does not evaluate shoulder expressions.

| wintype | Analysis a | Synthesis s |
|---|---|---|
| 0 | c | c |
| 1 | √c | a³ |
| 2 | 1 | c² |

Adjacent shoulders satisfy `aL*sL+aR*sR=1` in ideal arithmetic. 2D analysis multiplies vertical/horizontal windows; synthesis preserves horizontal-then-vertical order. Spatial inverse FFT already divides by block area. For O=1, cL=cR=√0.5, so two identity contributions of a unit sample are 0.5 each, summing to 1.

## DFTTest raw windows

At length L and index j, use half-sample coordinates `r=(j+0.5)/L`, `θ=2πr`, evaluated in binary64. Equations define behavior; names are labels. ID 9 has no abs(r−0.5): at L=2 it produces `[0.74,0.50]`, genuinely asymmetric. ID 4 uses a finite I0a approximation rather than an arbitrary Bessel implementation.

| ID | Formula |
|---|---|
| 0 | 0.5 − 0.5 cos θ |
| 1 | 0.53836 − 0.46164 cos θ |
| 2 | 0.42 − 0.5 cos θ + 0.08 cos 2θ |
| 3 | 0.35875 − 0.48829 cos θ + 0.14128 cos 2θ − 0.01168 cos 3θ |
| 4 | I0a(πβ√(1−(2r−1)²)) / I0a(πβ) |
| 5 | Σ(m=0..6) a[m] cos(mθ) |
| 6 | 0.2810639 − 0.5208972 cos θ + 0.1980399 cos 2θ |
| 7 | 1 |
| 8 | (2/L)(L/2 − abs(j+0.5−L/2)) |
| 9 | 0.62 − 0.48(r−0.5) − 0.38 cos θ |
| 10 | 0.355768 − 0.487396 cos θ + 0.144232 cos 2θ − 0.012604 cos 3θ |
| 11 | 0.3635819 − 0.4891775 cos θ + 0.1365995 cos 2θ − 0.0106411 cos 3θ |

```text
a = [0.27105140069342415, -0.433297939234486060,
     0.218122999543110620, -0.065925446388030898,
     0.010811742098372268, -7.7658482522509342e-4,
     1.3887217350903198e-5]
I0a(x):
    numerator = denominator = total = 1
    for k = 1..14:
        numerator *= x/2
        denominator *= k
        v = numerator/denominator
        total += v*v
        if v <= 1e-8: break
    return total
```


Negative finite coefficients are allowed. Active window overflow or zero normalization denominators fail creation.

## Forming the space/time volume

smode=0 keeps raw spatial sw. smode=1 uses step S=B−O and squared sums of equal-phase samples: `sw[j]=raw[j]/sqrt(sum raw[k]²)`. Summation visits j,j−S,..., then j+S,j+2S,... . tmode=0 keeps raw tw; tmode=1 performs analogous squared-sum normalization using the temporal step.

```text
V = T*B*B
h[z,y,x] = float(tw[z]*sw[y]*sw[x]/sqrt(V))
E = sum_float_in_z_y_x_order(h*h)
wscale = 1/E
G = FFT(255*h)
```

Products are formed in double before float storage; E accumulates stored float values sequentially. The mean template must be finite, with Re(G[0])≠0 when zmean=true. For example, flat-top B=4 with zero spatial overlap can normalize to a sign window with zero template DC. Disabling zmean permits that finite window; retaining zmean cannot divide by zero.

Analysis input is q*h, and U is the normalized inverse. tmode=0 precomputes `float(h*V)` and synthesizes `U*float(h*V)` at z=floor(T/2). tmode=1 first scales the inverse and synthesizes `(U*V)*h` for each covering block's target slice. The orders are algebraically equivalent but round differently in binary32. smode=0 emits the spatial center; smode=1 adds the whole block. Quantize once at the end, never per block.

B=2,O=1,T=1 with rectangular windows gives sw=1/√2,h=1/4,E=1/4,wscale=4. Each identity block contributes q/4, and four blocks reconstruct q. Type-0 sigma=8 gives threshold `8/wscale=2`, not 64.

For tmode=0, identity gain with spatial OLA is `tw[center]²`; center mode also multiplies by `sw[center]^4`. It is not universally 1: twin=6 has center value 1.000001. tmode=1 temporal phase normalization gives ideal total temporal gain 1 while retaining the spatial center-window gain. See [temporal overlap](../dfttest/temporal-ola.md).
