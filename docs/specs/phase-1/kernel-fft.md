# Real FFT adapter

Specification: FFT-001. Evidence: new adapter contract; normalization conversion is derived from both fixed CPU references.

## Inputs and outputs

Phase 1 implements positive 2D real shapes `[H,W]`, with independently transformable active blocks. Real storage is binary32; frequency storage consists of live `std::complex<float>` objects with logical shape `[H,K]`, K=floor(W/2)+1. Coordinates are y,x and ky,kx. Batch is an outer axis and never participates in a transform.

Descriptors carry explicit element strides/distances, capacity C and active count c. Require `0 <= c <= C`; c=0 is a no-access no-op and does not require dereferenceable buffer pointers. For c>0, all active slices meet [MEM-001](kernel-plane-geometry.md). Padding and inactive blocks are not inputs or outputs. No hidden extra real columns are required of callers.

Forward input and inverse input are logically read-only; input/output do not alias. If a backend destroys its C2R input or uses planning buffers destructively, the adapter supplies private mutable scratch. Plan creation cannot consume a caller's pixel data. Cache entries and precomputed grid spectra are immutable.

## Transform definition

Let N=H*W and i be the imaginary unit:

```text
X[ky,kx] = sum(y=0..H-1,x=0..W-1)
           a[y,x] * exp(-2*pi*i*(ky*y/H + kx*x/W))

a[y,x] = (1/N) * sum(ky=0..H-1,kx=0..W-1)
         X[ky,kx] * exp(+2*pi*i*(ky*y/H + kx*x/W))
```

Forward has scale 1; inverse has scale 1/N. N excludes batch count. The adapter applies inverse normalization exactly once. Both old projects use unnormalized inverse transforms, so each algorithm's synthesis contract explicitly accounts for this difference.

The missing half is `X[ky,kx] = conjugate(X[(-ky) mod H,(-kx) mod W])`. The kx=0 column, and kx=W/2 for even W, must themselves obey this Y-axis symmetry. They are not entirely real columns; only bins self-conjugate on **all** axes have zero imaginary component. Odd W has no Nyquist column. Do not manufacture a Nyquist bin or zero a whole boundary column.

Spectral kernels preserve Hermitian compatibility by symmetric real gains and compatible grid subtraction/restoration. The inverse's domain is a Hermitian-compatible spectrum; inconsistent arbitrary complex input is not a public filtering mode. Unit-test injection of such input must fail validation rather than define a second transform.

## Energy and frequency checks

For real a, Parseval is `sum(a*a) = sum(weight[kx]*abs(X)^2)/N` over the stored half. Weight is 1 for DC and even-width Nyquist columns, and 2 otherwise. This weighting is for total-energy verification. Neither algorithm doubles its per-bin noise threshold merely because the spectrum is compressed.

Full-axis signed bins are k for k <= floor(n/2), otherwise k-n; magnitude uses min(k,n-k). First-batch constant-noise filtering does not need a radial frequency approximation.

## Examples

- H=1,W=4, a=[1,0,0,0]: stored X=[1,1,1]; normalized inverse recovers a.
- H=1,W=4, a=[0,1,0,0]: stored X=[1,-i,-1]. This detects a reversed forward sign.
- H=1,W=3, a=[1,1,1]: stored X=[3,0], with no Nyquist. Energy is 9/3=3.
- H=2,W=3, all samples 2: X[0,0]=12, other bins zero; inverse samples are 2, not 12.
- Two identical blocks have identical spectra to one block; normalization must not acquire a factor of two. A final batch with c=1,C=4 leaves slots 1..3 untouched.

## Precision, failures and acceptance

PocketFFT is the required backend; FFTW is optional. Backend selection cannot change the logical contract. Use a fixed dependency revision, CPU ISA policy and one FFT thread per request in phase 1. No global fast-math or implicit FP contraction. Record backend-specific rounding in differential reports; bit identity across FFT implementations is not assumed.

Reject unrepresentable dimensions/distances and unsupported backend requests at creation. Reject non-finite input/output or transform overflow with a frame error. Backend plans/resources belong to the instance or exclusive workspace and follow [RUN-001](execution.md).

Verify forward sign, axis order, odd/even compression, non-square shapes, Parseval, normalized inverse, strided rows and partial batches against an independently written binary64 direct DFT. Round-trip alone is insufficient. The current specs do not authorize 3D execution; phase 2 extends the shape contract to `[T,H,W]` and defines N=T*H*W.
