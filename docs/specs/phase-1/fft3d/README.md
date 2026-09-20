# FFT3D specifications

Phase-1 spatial specifications are drafted below, organized into plugin interfaces and mathematical kernels. They are implementation contracts derived from the [reference version](../README.md#reference-versions), not runtime acceptance results.

| Specification | Phase-1 content |
| --- | --- |
| [plugin.md](plugin.md) | Complete VS argument order/defaults, admitted values, deferred-feature gates, output and errors |
| [kernel-geometry.md](kernel-geometry.md) | Per-plane block grid, reflected cover, safe domain and visible crop |
| [kernel-window-reconstruction.md](kernel-window-reconstruction.md) | Three windows, pixel/noise units, overlap reconstruction and integer/float output |
| [kernel-wiener.md](kernel-wiener.md) | Constant-noise spatial Wiener, beta floor, degrid template and restoration |

Phase 1 requires bt=1 and includes default degrid=1, all three windows and scalar/Highway VS comparison. Retain public bt default 3 and explicitly reject it until temporal processing is implemented. Empty/omitted planes selects all; float output clipping includes chroma. See [acceptance](../acceptance.md).
