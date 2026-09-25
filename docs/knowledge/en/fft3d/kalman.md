# FFT3D Kalman state and random access

[Contents](../README.md)

bt=0 recurs over every spatial block/bin instead of using a finite temporal FFT window. Output 0 bypasses Kalman and copies source data except for explicitly uninitialized AVS mode=1 planes. Virtual S0 is initialized directly; source frame 0 is not transformed to seed it. S_n means state after consuming source n.

[FFT3D API](../../../api/en/fft3d.md) · [Noise and enhancement](noise-enhancement.md)

## One recurrence step

L is a complex estimate, C covariance and Q process covariance. X is the current centered/windowed spatial FFT, without recurrence degrid subtraction. Uniform noise uses R=R0; analytic/sampled noise uses R=max(P,1e−15). Positive pfactor's magnitude does not scale the Kalman sampled table; sigma still initializes C/Q.

```text
R0 = sigma_eff^2 / (1/(bw*bh))
S0: L=0+0i, C=Q=R0
motion = (X.re-L.re)^2 > R*kratio^2
      or (X.im-L.im)^2 > R*kratio^2
if uniform R==0:
    L=X; C=Q=0
elif motion:
    L=X; C=Q=R
else:
    sum = C+Q
    gain = sum/(sum+R)
    Q_new = (gain*gain)*R
    C_new = (1-gain)*sum
    L_new = gain*X+(1-gain)*L
```


Motion compares squared real/imaginary differences separately; either exceeding the threshold resets both components. Equality smooths. This is not squared complex magnitude. Preserve the stated binary32 operation order; enhanced spectra, inverse output and quantized pixels never feed back into state.

## Two-step example

R0=R=2,kratio=2,X1=1+0i gives threshold 8 and squared difference 1, so sum=4,gain=2/3,L1=2/3,C1=4/3,Q1=8/9.

For X2=1+0i, motion still does not trigger. sum=20/9,gain=10/19, giving L2=16/19,C2=20/19,Q2=200/361. Changing the first input to X1=3+0i gives 9>8 and resets to L1=3,C1=Q1=2. These fractions illustrate ideal arithmetic; do not replace the prescribed float steps with simplified expressions.

Logical real/imaginary C components remain equal, as do Q components. Each bin therefore stores complex L and two distinct float C/Q values: 16 bytes rather than six floats/24 bytes. C and Q are not equal and cannot be merged with each other.

## Restoring random requests

With W=kalman_warmup (default 8), define earliest=max(0,n-W-1). At n>0 request start, freeze the largest completed checkpoint j in [earliest,n]. Restore a private copy, or initialize L=0,C=Q=R0 at virtual j=earliest if no eligible checkpoint exists. Consume j+1..n in order, at most W+1 recurrence frames. An exact checkpoint j=n only needs source n for properties/copied regions. Frame 0 does not advance or reset state.

At n=10000,W=8: cold consumes 9992..10000; checkpoint 9996 needs 9997..10000; checkpoint 9991 is eligible at the budget edge, but 9990 is too old. W=0 permits only the current recurrence frame on a cold request, yet a checkpoint at n-1 still preserves continuous long history. A backward request cannot use a future checkpoint. No implicit full-prefix fallback occurs after eviction.

Rendering copies L_n into a private spectrum, then enhances, inverts, synthesizes and clips. ROI/field packing belongs to instance configuration; instances do not share state.

## Budget and cost

The default cache budget is `max(64 MiB, actual capacity of one complete checkpoint plus metadata)`, with at most four checkpoints and LRU eviction. A state exceeding 64 MiB can still retain one checkpoint. This is not a fixed 64 MiB hard limit and is not controlled by cache_mb/cache_frames.

For 1280×720 YUV420, 32×32 blocks, overlap 16 and three planes, state payload is about 47.265 MiB, versus 70.897 MiB for expanded six-component state. Active private states, starting snapshots, scratch and host caches consume additional memory. Leased evicted snapshots can remain alive.

For fixed W, cold recurrence work is bounded independently of n. A sampled model may first fetch pframe outside the window and can reacquire that frame during recurrence; W+1 is not a total callback or wall-clock bound. Staged release limits retained source ownership; upstream caches and host metadata remain outside plugin memory guarantees. Failed work publishes no checkpoint.

## History dependence and concurrent requests

Only completed output states are published; intermediate warmup states are not cached. Requests privately advance immutable leased snapshots and never wait for another output. Concurrent cold requests may duplicate work and select different histories. The first published state for a frame wins while retained; another in-flight request still returns its own result. After eviction, recomputing the frame can choose another history and differ. This is a deliberate change from unlimited canonical replay, not a claim of request-order independence.

Continuous processing normally extends the previous state one frame at a time. Warmup is a work budget, not a rolling history cap. More warmup can improve startup but does not guarantee rapid recovery of long-history denoising. Motion resets operate per frequency bin, not as scene-boundary detection. See the [public parameter](../../../api/en/fft3d.md).
