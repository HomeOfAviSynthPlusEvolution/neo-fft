# DFTTest specifications

Phase-1 spatial specifications are drafted below, organized into plugin interfaces and mathematical kernels. The [reference version](../README.md#reference-versions) establishes static evidence; black-box results remain pending.

| Specification | Phase-1 content |
| --- | --- |
| [plugin.md](plugin.md) | Full VS argument order/defaults, normalization, admitted values and deferred-feature errors |
| [kernel-geometry.md](kernel-geometry.md) | Center-sample/overlap-add grids, reflection and crop |
| [kernel-window-reconstruction.md](kernel-window-reconstruction.md) | All twelve windows, energy calibration, input scale, synthesis and no-dither conversion |
| [kernel-filter.md](kernel-filter.md) | zmean, ftype 0..4, f0beta branches, thresholds and gains |

Phase 1 requires tbsize=1 and includes both spatial modes, all five ftypes and all windows, with zmean on/off and no dither. Public tbsize still defaults to 3 and is explicitly rejected until supported. Sigma2 independently defaults to 8; empty planes selects none; window 9 follows the asymmetric source equation. See [acceptance](../acceptance.md).
