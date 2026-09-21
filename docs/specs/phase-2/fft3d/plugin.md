# FFT3D: phase-2 VS interface addendum

Specification: F3D-VS-002. Keep [F3D-VS-001](../../phase-1/fft3d/plugin.md) registration, parameter order, types, defaults and restrictions, with the following overrides:

- bt defaults to 3 and now accepts every integer 1..5. Reject bt=-1,0 as deferred enhancement/Kalman modes and all values outside 1..5. T=2 and 4 are valid; no odd-length restriction applies.
- Default bt=3 is now usable. A short clip is legal: each request whose full neighborhood is unavailable uses effective btcur=1, without changing the configured bt or a different request's plan.
- degrid now follows [F3D-TIME-001](kernel-temporal.md) for btcur>1: use current-frame block DC and subtract T times its grid only from temporal frequency zero. The old phase-1 per-block rule remains the btcur=1 case.

sigma2..4 must still equal parsed sigma; pattern controls retain phase-1 default restrictions; sharpen=dehalo=0, interlaced=false, zero ROI, mt=false and existing opt/ncpu/backend policy remain. No implicit enablement of phase-3/4 modes is allowed. Sigma/beta/degrid domains and sample-format scaling are unchanged.

At creation validate selected-plane geometry, immutable spatial windows/grid, finite noise constants for possible effective sizes, checked sizes and execution controls. At evaluation select btcur and slots by [RUN-TIME-001](../execution.md), transform raw spatial blocks, apply temporal Wiener/degrid, reconstruct cur and run the existing spatial synthesis/conversion. Never filter neighbor cached spectra in place.

Source properties and copied unselected planes come from clip[n], never the first/last neighbor. Output geometry, format, frame count/rate are unchanged. Omitted or empty planes still selects all, unlike DFTTest. Runtime failures follow phase-1 atomic output, owned error and cleanup contracts.

Examples: N=1,bt=5 is legal and all output uses btcur=1; N=5,bt=4 uses btcur=4 at n=2,3 and fallback at n=0,1,4. A fallback frame must match the explicit bt=1 call within the same implementation's deterministic execution contract. This intentional boundary fallback does not authorize silently replacing an invalid configured bt with 1.
