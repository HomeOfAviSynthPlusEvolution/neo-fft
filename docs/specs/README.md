# Implementation specifications

[Phase 1](phase-1/README.md): Shared FFT infrastructure, spatial FFT3D (`bt=1`), and spatial DFTTest (`tbsize=1`).

[Phase 2](phase-2/README.md): Stateless temporal filtering, multi-frame FFT3D (`bt=2..5`), and temporal DFTTest (odd `tbsize=3..15`). Includes public interface addenda, 3D FFT layout, window-energy/PSD/mean contracts, scheduling and acceptance.

[Phase 3](phase-3/README.md): Noise models and spectral enhancement: FFT3D frequency-dependent sigma, sampled patterns, preview, sharpen/dehalo and `bt=-1`; DFTTest curves and sampled noise. Includes precedence, calibration, dependencies, explicit reference differences and acceptance.

[Phase 4](phase-4/README.md): Complete VS CPU behavior: deterministic FFT3D Kalman replay with bounded checkpoints, ROI and field packing, DFTTest dither, execution controls and resource lifetime. Includes integration of the fixed DS2 staged-frame release API and acceptance criteria.

[Phase 5](phase-5/README.md): Original DFTTest temporal block overlap-add (`tmode=1`), even/odd temporal sizes, absolute block lattice, calibrated windows/profiles, deterministic execution and comparison with the pinned original AVS implementation. AVS integration and final distribution move to phase 6.

Current execution policy, revised through `c761c8c`, is defined in [RUN-004](phase-4/execution.md): both filters execute inline on the calling host thread, while host frame concurrency remains supported. DFTTest threads is reserved and has no execution or retention effect. FFT3D mt/ncpu/measure and DFTTest fft_threads are removed; KernelInfo retains its read-only fft_threads=1 diagnostic. Earlier phase signatures and admission rules are historical baselines, superseded by these explicit revisions.

Each function directory contains its mathematical operators in `kernel-*.md` and its public interface in `plugin.md`. Shared operators are referenced rather than duplicated. Implement the scalar definitions before adding optimized paths.
