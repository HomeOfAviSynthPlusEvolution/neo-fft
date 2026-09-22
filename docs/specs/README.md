# Implementation specifications

[Phase 1](phase-1/README.md): Shared FFT infrastructure, spatial FFT3D (`bt=1`), and spatial DFTTest (`tbsize=1`).

[Phase 2](phase-2/README.md): Stateless temporal filtering, multi-frame FFT3D (`bt=2..5`), and temporal DFTTest (odd `tbsize=3..15`). Includes public interface addenda, 3D FFT layout, window-energy/PSD/mean contracts, scheduling and acceptance.

[Phase 3](phase-3/README.md): Noise models and spectral enhancement: FFT3D frequency-dependent sigma, sampled patterns, preview, sharpen/dehalo and `bt=-1`; DFTTest curves and sampled noise. Includes precedence, calibration, dependencies, explicit reference differences and acceptance.

Each function directory contains its mathematical operators in `kernel-*.md` and its public interface in `plugin.md`. Shared operators are referenced rather than duplicated. Implement the scalar definitions before adding optimized paths.
