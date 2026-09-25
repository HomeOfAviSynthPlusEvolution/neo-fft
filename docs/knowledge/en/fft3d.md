# FFT3D: from pixels to spectra and back

[Contents](README.md)

FFT3D transforms overlapping image blocks, suppresses frequencies according to estimated noise, and combines the reconstructed blocks. bt=2..5 also transforms across frames at the same positions; there is no motion alignment. See the [complete API](../../api/en/fft3d.md).

## 1. Choose the processing path

```mermaid
flowchart LR
    A[Selected plane ROI] --> B[Optional field packing and reflection]
    B --> C[Center samples and apply analysis window]
    C --> D[Spatial FFT]
    D --> E[Optional temporal FFT]
    E --> F[Grid compensation and Wiener]
    F --> G[Recover current spatial spectrum]
    G --> H[Mode-specific enhancement]
    H --> I[Spatial inverse FFT and overlap-add]
    I --> J[Write ROI and quantize]
```

Uniform single-frame processing fuses Wiener and enhancement using the same prefilter residual power. The diagram describes functional relationships, not a universal sequence of separate enhancement calls. Effective pshow directly reconstructs a selected preview block; bt=-1 skips denoising; bt=0 substitutes Kalman recurrence for temporal FFT/Wiener. Parameters select modes at creation; temporal endpoint fallback depends on output position.

## 2. Blocks and the working image

W,H denote the current plane's ROI size, Bx,By the block size. Chroma uses the same block dimensions in its own samples. On either axis:

$$S=B-O,\qquad N=\left\lceil\frac{W-O}{S}\right\rceil+2,\qquad P=NS+O.$$

The ROI occupies [S,S+W) in the padded image, and block i starts at iS. Two extra blocks supply full edge contributions. For padded coordinate x, j=x−S maps to:

$$r(j)=\begin{cases}-j&j<0\\j&0\le j<W\\2W-2-j&j\ge W.\end{cases}$$

This is one reflection without duplicated endpoints. Creation requires each padding extent to be less than the corresponding ROI dimension; reflection is not repeated indefinitely.

For a 64×48 ROI, 8×6 blocks and 2×2 overlap, steps are 6×4, the grid is 13×14 and the padded image is 80×58 with ROI origin (6,4). Reconstruction crops back to the same ROI, preserving output dimensions.

interlaced=True permutes height-8 ROI rows to [0,2,4,6,7,5,3,1], filters that full-height image, then undoes the permutation. It does not split independent half-height fields or change temporal frame numbering.

## 3. Sample, window and spectrum scales

Integer samples remain in their native code-value domain. Only integer YUV chroma subtracts 2^(bits−1); luma, RGB and all float planes have no midpoint subtraction. Neutral 10-bit chroma 512 becomes zero, whereas RGB 512 stays 512. Public sigma is an 8-bit standard deviation:

$$\sigma_e=\begin{cases}\sigma\,2^{bits-8}&\text{integer}\\\sigma/255&\text{float32}.\end{cases}$$

sigma=2 gives 8 at 10-bit and about 0.00784314 for float. Float input pixels themselves are not multiplied by 255.

Nonoverlap window regions equal 1. For overlap shoulder j=0..O−1:

$$c_L(j)=\cos\frac{\pi(j-O+1/2)}{2O},\qquad c_R(j)=\cos\frac{\pi(j+1/2)}{2O}.$$

| wintype | Analysis a | Synthesis s |
|---:|---|---|
| 0 | c | c |
| 1 | √c | a³ |
| 2 | 1 | c² |

Analysis successively computes (sample−base)·ay·ax. Synthesis multiplies the inverse transform by its window. In ideal arithmetic aL·sL+aR·sR=1, so overlapping contributions add without per-pixel weight division. These paired windows and their binary32 tables cannot be replaced with arbitrary named windows. O=0 uses all ones without evaluating shoulder formulas.

The spatial forward FFT uses the negative exponent and no normalization; real input retains By×(floor(Bx/2)+1) complex bins. Inverse normalization is 1/(Bx·By). Reconstruction does not divide by area again or double half-spectrum power.

## 4. Wiener retention

Let X be a spatial block spectrum and G the spectrum of a full-scale constant block multiplied by the analysis window. Grid compensation computes:

$$g=degrid\,\frac{\operatorname{Re}X_{0,0}}{\operatorname{Re}G_{0,0}},\qquad M_k=gG_k,\qquad R_k=X_k-M_k.$$

degrid=0 gives M=0,R=X. Because G contains window-induced frequencies, this removes more than DC. At degrid=1 an ideal constant block belongs entirely to M; restoring it preserves brightness.

Uniform spatial noise power uses block area, not actual window energy:

$$\nu=B_xB_y\sigma_e^2,\qquad P=\operatorname{Re}(R)^2+\operatorname{Im}(R)^2+10^{-15},$$
$$w=\max\left(\frac{P-\nu}{P},\frac{\beta-1}{\beta}\right),\qquad Y=wR+M.$$

beta=1 suppresses residual power at or below noise to zero; beta=2 retains at least half the residual amplitude. This compares powers, not complex amplitude against sigma.

Model priority is positive pfactor sampling, then unequal sigma1..4 profile, then uniform noise. The analytic profile interpolates sigma before squaring and spatial FFT calibration. Sampling instead directly estimates residual power:

$$P_k^{sample}=pfactor\,|X_k-gG_k|^2\,\frac{f_x^2+f_y^2}{f_x^2+f_y^2+pcutoff^2},$$

where fx=2x/Bx and fy=2min(y,By−y)/By. The analytic profile uses different coordinates and interpolation nodes; see [noise and enhancement](fft3d/noise-enhancement.md). Automatic sampling searches bx=2..Nx−3, by=2..Ny−3, requiring at least 5×5 blocks. It minimizes weighted power before pfactor, selecting the first minimum in row-major order; manual indices directly select the block. Sampled power must not be squared again or multiplied by block area again.

## 5. Multiple frames, one current output

T is the request's effective temporal length. An incomplete configured neighborhood falls back wholly to T=1.

| bt | Source offsets, earliest first | Current temporal index c |
|---:|---|---:|
| 2 | −1, 0 | 1 |
| 3 | −1, 0, +1 | 1 |
| 4 | −2, −1, 0, +1 | 2 |
| 5 | −2, −1, 0, +1, +2 | 2 |

For one spatial bin, transform T values without another temporal analysis window:

$$F_{m,k}=\sum_{j=0}^{T-1}X_{j,k}\exp(-2\pi i jm/T).$$

Uniform power becomes T·Bx·By·sigma_e², and profile/sampled tables are multiplied by T. degrid takes the **current block's** DC ratio and subtracts T·M only from temporal frequency m=0. After per-frequency Wiener and grid restoration, recover only current index c:

$$Y_k=\frac1T\sum_{m=0}^{T-1}F^{out}_{m,k}\exp(+2\pi i cm/T).$$

Then perform spatial inverse FFT and synthesis. For T=2 in the stated order, current output is (Fout0−Fout1)/2, not their sum. Adjacent outputs reuse some X, but current DC, neighborhoods and filtered Y differ. The [spectrum cache](fft3d/spectra-cache.md) stores raw X.

## 6. Kalman and enhancement placement

Kalman skips the temporal DFT. Each spatial bin retains estimate L, covariance C and process covariance Q. Virtual initialization is L=0,C=Q=R0, where R0 uses sigma and block area. Source frame 0 does not initialize state, and output frame 0 is copied.

Starting with source frame 1, if either squared real or imaginary change exceeds R·kratio², reset both parts to current X and C,Q to R. Otherwise, using old C,Q:

$$s=C+Q,\quad k=\frac{s}{s+R},\quad L'=kX+(1-k)L,\quad C'=(1-k)s,\quad Q'=k^2R.$$

Uniform R=0 has an identity branch avoiding 0/0. Profile/sampled models use positive per-bin floors; positive Kalman pfactor selects the model without scaling its power. Random output n uses the nearest completed checkpoint whose continuation costs at most kalman_warmup+1 steps, or initializes at max(1,n-kalman_warmup). The default is 8 preceding frames plus the target; sequential requests extend long history. Cache/request history and concurrency can change results.

Real/imaginary covariance components remain equal within C and within Q. Storage therefore needs one float C, one float Q and complex L: 16 bytes per bin. C and Q are distinct quantities. See [Kalman state and replay](fft3d/kalman.md).

Enhancement uses spatial residual power q and frequency weights. Its position depends on the path:

| Path | Spectrum/power used for enhancement |
|---|---|
| T=1, uniform noise | Fused with Wiener, sharing original residual power |
| T=1, profile/sampled | Wiener first, then recompute enhancement residual and power |
| T>1 | After temporal Wiener and current spatial-spectrum recovery |
| bt=-1 | Original spatial spectrum |
| Kalman | Private copy of L; enhancement never feeds back into recurrence |

Zero sharpening and dehalo strengths skip enhancement arithmetic. Effective pshow overrides all denoising/enhancement.

## 7. Complete 4×4 example

Take a 16×16 GRAY8 image whose rows repeat [10,14]. Use bw=bh=4, ow=oh=0, bt=1, sigma=4, beta=1, degrid=0 and no enhancement. The image satisfies padding geometry and the window is rectangular. Each block is:

```text
10 14 10 14
10 14 10 14
10 14 10 14
10 14 10 14
```

The unnormalized 2D FFT has only DC=192 and horizontal Nyquist=−32. Noise power is 16×4²=256. Ignoring the negligible epsilon:

| Bin | Input | Power | Gain | Output |
|---|---:|---:|---:|---:|
| DC | 192 | 36864 | 143/144 | 190⅔ |
| Horizontal Nyquist | −32 | 1024 | 3/4 | −24 |

Inverse FFT divides by 16:

$$u_{even}=\frac{190\tfrac23-24}{16}=10.416\overline6,\qquad u_{odd}=\frac{190\tfrac23+24}{16}=13.416\overline6.$$

Without overlap each pixel has one contribution. Add 0.5, clip and truncate to obtain repeated **[10,13]**. The difference shrinks from 4 to 3 and the mean drops slightly because degrid was deliberately disabled, exposing DC to Wiener. Default degrid=1 removes/restores the constant component and does not follow this direct DC attenuation example.

## 8. Reconstruction, boundaries and determinism

The spatial inverse is already normalized. Blocks receive horizontal synthesis weights and accumulate left to right, then vertical weights and accumulation top to bottom. This canonical order also applies under host concurrency; arbitrary floating-point atomic accumulation would change results.

Final integer conversion forms (reconstructed+0.5)+base, clips to the bit-depth range and safely truncates toward zero. Float output clips directly to [0,1], **including float chroma**. Copies of unselected planes, areas outside the ROI and Kalman frame 0 skip quantization. Explicit AVS mode 1 skips output writes entirely, including frame 0.

FFT3D creates no workers. Selected planes and Kalman replay execute in the calling thread; host requests may run concurrently with separate workspaces and mutable replay state. Within one arithmetic configuration, caching, request order and scheduling should not change non-Kalman output. Bounded Kalman warmup does not provide that guarantee. Changing scalar/Highway or FFT paths can change rounding, especially near gain thresholds, reset thresholds and half-LSB boundaries. One bitwise-equal example is not a universal equality guarantee.

Unused arguments still undergo their specified domain checks. Consumed inputs and active intermediates must be finite; geometry, allocation or numerical errors fail creation/frame processing instead of silently substituting algorithms. See [sample scales](shared/sample-domain.md), [windows](shared/windows-reconstruction.md) and [execution precision](shared/execution-precision.md).
