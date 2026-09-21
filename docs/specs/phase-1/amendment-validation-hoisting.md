# Phase-1 Amendment: Validation Hoisting and Precondition Guarantees

Specification: AMEND-001. Amends: [RUN-001](execution.md) and [FFT-001](kernel-fft.md). Status: Approved specification amendment.

## 1. Background and Purpose

The initial Phase-1 specification located defensive checks—including non-finite sample scanning, memory boundary verification, and boundary spectral symmetry assertions—inside the lowest-level spectral transform operators.

Because spatial filtering decomposes an image plane into thousands of overlapping blocks per frame, executing redundant sample-by-sample floating-point scans and boundary symmetry assertions inside each individual block transform introduces an unnecessary computational bottleneck.

This amendment formally defines a two-tier contract structure: external untrusted inputs are validated strictly at the pipeline admission boundary, while inner computational operators execute under guaranteed mathematical and memory preconditions.

## 2. Invariant and Boundary Realignment

### 2.1 Pipeline Admission Boundary (Untrusted to Trusted Boundary)
All verification of external, untrusted host data occurs at the frame evaluation boundary before spatial decomposition begins:
- **Geometry and Memory ([MEM-001](kernel-plane-geometry.md))**: Plane dimensions, row strides, memory alignment, and non-overlapping buffer extents are verified once per selected plane.
- **Finite Value Verification ([RUN-001](execution.md))**:
  - For `Float32` selected planes, all visible input samples must be verified finite (`NaN` and `+/-Inf` rejected) before block gathering commences. Encountering any non-finite sample immediately terminates frame evaluation with a descriptive `'non-finite'` frame error.
  - For integer planar formats (`UInt8`, `UInt10`, `UInt12`, `UInt14`, `UInt16`), visible sample values are bounded integers by definition; floating-point finite scanning is not required.
- **Observable Error Behavior**: External observable error behavior remains strictly identical to [RUN-001](execution.md). Invalid parameters and non-finite inputs continue to fail deterministically, and unselected planes continue to copy bitwise.

### 2.2 Inner Operator Precondition Contract (Trusted Domain)
Inner computational operators—including window application, discrete spectral transforms, frequency-domain gain modulation, and spatial overlap reconstruction—operate exclusively within memory buffers managed by the pipeline.
- **Preconditions Assumed**:
  - Buffer pointers are valid, correctly aligned, and disjoint.
  - Input sample values are guaranteed finite by the admission boundary.
  - Forward real-to-complex transforms produce a valid half-spectrum obeying Hermitian symmetry across boundary columns (DC and Nyquist).
  - Frequency-domain filtering modulates spectral bins using real-valued gains, preserving boundary conjugate symmetry by mathematical invariant.
- **Elimination of Inner Per-Block Scans**:
  - Lowest-level spectral transform operators shall not perform redundant per-block sample scanning, buffer re-validation, or boundary symmetry tests.
  - Inner operators execute mathematical transforms directly, assuming all preconditions are satisfied by the calling pipeline.

## 3. Uniform Behavior

- The reallocation of validation responsibility applies unconditionally across all compilation targets and build configurations.
- There is no distinct diagnostic vs. release divergence or configurable runtime toggle for inner operator validation. The contract structure is unified: admission verifies external data, and inner kernels execute under verified preconditions.

## 4. Acceptance and Verification Impact

- **Black-box Acceptance ([P1-ACCEPT-001](acceptance.md))**: Admission tests verifying the rejection of non-finite inputs on selected planes must pass without modification.
- **Numerical Invariance**: The hoisting of verification does not alter numerical calculation, arithmetic precision, or rounding; differential comparison tolerances remain unchanged.
- **Operator Verification**: Standalone tests of the spectral transform operator shall supply inputs satisfying its defined precondition domain.
