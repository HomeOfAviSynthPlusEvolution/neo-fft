# FFT3D: phase-3 VS interface

Specification: F3D-VS-003. Extends [phase-2 interface](../../phase-2/fft3d/plugin.md). Defaults/types remain those of phase 1 unless stated here.

| Parameter | Default | Phase-3 domain/effect |
| --- | --- | --- |
| bt | 3 | -1 or 1..5; -1 is spatial enhancement only; 0 remains unsupported |
| sigma2, sigma3, sigma4 | Each inherits parsed sigma | Finite >=0; frequency-dependent noise |
| pfactor | 0 | Finite >=0; positive selects sampled pattern and scales power linearly |
| pframe | 0 | int32, clamped to [0,N-1] after N validation |
| px, py | 0,0 | Block-grid indices, not source pixels; both zero means automatic search |
| pcutoff | 0.1 | Finite >0; sample high-pass weight |
| pshow | false | Preview selected noise block when effective pattern mode is enabled |
| sharpen | 0 | Finite >=0, no upper clamp |
| scutoff | 0.3 | Finite >0 |
| svr | 1 | Finite >=0, including zero |
| smin, smax | 4,20 | Finite 0<=smin<=smax; zero maximum yields no sharpen gain |
| dehalo | 0 | Finite >=0, no upper clamp |
| hr | 2 | Finite >0; active dehalo requires a finite positive sampled window maximum |
| ht | 50 | Finite >=0; historical unscaled amplitude, see enhancement contract |

Other phase-2 restrictions remain, including kratio=2, interlaced=false, zero ROI, mt=false and opt/backend controls. These nonnegative/positive domains are the project's safe admitted domain, not claims that old code validates them. Validate input domains even if inactive. Validate px/py>=0 always; for active sampling/preview, bounds must fit every selected plane's grid. Automatic search requires Nx>=5, Ny>=5. Dormant coordinates need not fit a grid.

## Selection and precedence

Resolve values in binary32. When pfactor=0, include existing format scaling of sigma1..4 before testing equality and require the scaled values to be finite. Let varying mean a scaled sigma2..4 differs from scaled sigma. Positive pfactor selects sampled mode directly; unused scaled sigmas need not be constructed. Consumed versus inactive derived tables follow [execution](../execution.md).

| Condition, in order | Effective model |
| --- | --- |
| pfactor>0 | Sampled pattern; sigma1..4 do not contribute to Wiener noise |
| pfactor=0 and varying | Analytic sigma profile, effective pattern strength 1 |
| pfactor=0 and all equal | Uniform sigma, pattern disabled |

If pshow=true and effective pattern mode is enabled, preview overrides bt and enhancement, selecting a block from n, not pframe. This also applies to analytic sigma mode. Otherwise pshow is inactive. bt=-1 without effective preview skips denoising and noise estimation regardless of noise controls. This dependency pruning is intentional; validate active preview geometry independently of bt.

Normal sampled denoising samples pframe per selected plane, using that plane's format, geometry and analysis window, and publishes immutable profiles. Analytic profiles need no additional frames. Noise controls do not select planes or change temporal geometry.

Output dimensions/rate/count/properties come from source n. Unselected planes are copied exactly. Preview follows [noise contract](kernel-noise.md), without text overlays or new frame properties.
