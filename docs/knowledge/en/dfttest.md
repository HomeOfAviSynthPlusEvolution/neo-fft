# DFTTest: space/time spectra and reconstruction

[Contents](README.md)

DFTTest forms a volume from T colocated B×B patches, filters residual spectra by power, restores the mean and reconstructs output. Defaults are B=16, spatial overlap 12, T=3 and tmode=0. There is no motion compensation or cross-frame raw-spectrum cache.

[Complete API](../../api/en/dfttest.md) · [Spectral models](dfttest/spectral-models.md) · [Temporal overlap](dfttest/temporal-ola.md)

## Input to output

1. Determine logical frame slots from temporal mode and request unique source frames. Active nlocation adds its sample intervals.
2. Reflect spatial boundaries per actual plane and convert samples to 8-bit amplitude units.
3. Multiply by space/time window h, transform volume V=T×B² and store half-spectrum `[T,B,floor(B/2)+1]`.
4. With zmean, compute `g=Re(X[0])/Re(G[0])`, G=FFT(255h), and residual R=X−gG.
5. Build primary parameter A from the model, apply ftype to R, then restore gG.
6. Apply normalized inverse FFT and select the target temporal slice. Synthesize `U*float(h*V)` with a precomputed window in tmode=0, or `(U*V)*h` in tmode=1, emitting spatial centers or adding whole overlapping blocks.
7. After all temporal contributions, crop and convert or dither once. Properties come from n.

The mean is derived from the entire volume's DC with window weighting, not independently from each frame. FFT epsilon and window-energy calibration depend on analysis scaling; moving 1/√V requires more than retaining the old thresholds.

## Two temporal meanings

tmode=0 requires odd T and slots `clamp(n−floor(T/2)+z,0,N−1)`; only center z is reconstructed. tmode=1 uses a fixed start lattice, permits even T and adds the slice at n from every covering block. It is not repeated invocation of centered mode.

Spatial smode=0 gathers an odd block around every pixel and emits its center. smode=1 uses step B−O and squared-overlap-normalized windows. See [block geometry](shared/block-geometry.md).

## Reproducible mean example

Choose B=1,smode=0, rectangular spatial/temporal windows, T=3,zmean=true,ftype=2,sigma=0. Three single-pixel amplitudes in 8-bit units are `[1,2,6]`.

V=3,h=1/√3. Since the window is constant, the DC-derived mean is `(1+2+6)/3=3`. Type-2 zero gain clears residuals. Restoring the mean, applying normalized inverse and Vh synthesis yields center output 3, although the original center was 2. With zmean=false the output is 0.

Float input `[1/255,2/255,6/255]` yields center `3/255`, without clipping. This checks the temporal-volume mean, type-2 gain semantics and inverse scale compensation together.

## Where controls enter

| Controls | Computation affected |
|---|---|
| sbsize/smode/sosize, tbsize/tmode/tosize | Origins, dependencies, target slices and contribution counts |
| swin/twin, sbeta/tbeta | Analysis/synthesis windows and power calibration |
| sigma, curves, nlocation/alpha | Primary thresholds/gain table; sampling overrides curves |
| ftype/f0beta, sigma2/pmin/pmax | Per-bin gain function |
| zmean | Removal/restoration of the window-shaped mean |
| dither/seed | Final UInt8 quantization, not spectra |

Types 0/1 use a power table; 2/3/4 use gain parameters, unlike FFT3D's squared standard deviation. Sampling affects 0/1 only; curves can affect all five types.

## Boundaries and precision

Temporal clamping keeps distinct logical window positions for duplicate frames. Sample rectangles must fit the source plane; output padding cannot extend them. planes=[] still validates configuration before copying n. Active samples must be finite; integers clip at final output, while float retains negative chroma and out-of-range values.

Fixed FFT groups may batch eight blocks; batching is not threading. All blocks currently run on the calling thread, synthesizing by temporal start, then Y, then X. Within one arithmetic configuration, request order, host concurrency and compatibility argument values such as threads must not change results. Cross-SIMD/FFT comparisons use scoped numerical budgets; see [execution precision](shared/execution-precision.md).
