# neo-fft knowledge base

[Contents](README.md)

Follow pixels through spectra, filtering and reconstruction in the current implementation. Full signatures, defaults and calling errors belong to the [API reference](../../api/en/README.md). Current host entries are VapourSynth and AviSynth. Both filters are internally single-threaded; the host may request frames concurrently.

## Start with a filter

| Article | Scope |
|---|---|
| [FFT3D](fft3d.md) | Spatial/temporal Wiener, ROI, field packing, Kalman and enhancement; complete 4×4 example |
| [DFTTest](dfttest.md) | Space-time volumes, mean removal, gains, center output and overlap-add; worked examples |
| [KernelInfo](kernel-info.md) | Automatic dispatch, FFT profiles, lanes versus threads and diagnostic scope |

## Shared foundations

| Article | Questions answered |
|---|---|
| [Sample domains](shared/sample-domain.md) | Bit depth, float, chroma midpoint, sigma scales and quantization |
| [Block geometry](shared/block-geometry.md) | Each algorithm's grids, padding, reflection and visible regions |
| [FFT conventions](shared/fft-conventions.md) | Signs, half-spectrum layout, DC/Nyquist, power and normalization |
| [Windows and reconstruction](shared/windows-reconstruction.md) | Paired FFT3D shoulders, 12 DFTTest windows, energy and overlap-add |
| [Execution and precision](shared/execution-precision.md) | Host concurrency, workspaces, accumulation order, finiteness and comparison scope |

## FFT3D topics

| Article | Scope |
|---|---|
| [Noise and enhancement](fft3d/noise-enhancement.md) | Uniform/profile/sample priority, degrid, Wiener, sharpening, dehalo and preview |
| [Kalman state and replay](fft3d/kalman.md) | Two-step recurrence, initialization, resets, checkpoints, random access and state memory |
| [Row spectrum cache](fft3d/spectra-cache.md) | Raw spectra, 128 MiB default, demand-based eviction, sharing and 1080p/4K examples |

## DFTTest topics

| Article | Scope |
|---|---|
| [Spectral models](dfttest/spectral-models.md) | ftype=0..4, curves/sampling, sigma calibration, zmean and threshold examples |
| [Temporal overlap-add](dfttest/temporal-ola.md) | Absolute time grid, even blocks, endpoint clamping, window normalization and accumulation |
| [Dither](dfttest/dither.md) | Integer quantization, error diffusion, random perturbation, seed and request-order determinism |

Users can begin with the filter overview and follow relevant topics. Maintainers can read sample domains → FFT conventions → windows/reconstruction → block geometry → filter overviews → execution/precision. Formula examples distinguish ideal arithmetic from binary32; cache miss ratios do not directly predict total speedup.

There is currently no DFTTest cross-frame spectrum cache or selectable FFT backend. Legacy execution arguments such as threads are accepted only by the host and completely ignored by the plugin. See [API migration](../../api/en/README.md#threads-and-execution-parameters) for compatibility arguments and the [migration guide](../../api/en/migration.md) for required script edits and differences in results.
