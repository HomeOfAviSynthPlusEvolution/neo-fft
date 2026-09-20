# Implementation specifications

[Phase 1](phase-1/README.md): shared FFT infrastructure, spatial FFT3D and spatial DFTTest, including VS differential acceptance.

Status: first-batch specification draft; no filter implementation or VS acceptance results yet.

Organize specifications by phase and function. Each function directory contains its mathematical operators in `kernel-*.md` and its public interface in `plugin.md`. Common contracts and shared operators belong to the phase that introduces them and are referenced rather than duplicated. Implement the scalar definitions before adding optimized paths.

Implementation and review requirements are defined within this specification tree. Reference versions are listed in the [phase-1 introduction](phase-1/README.md#reference-versions). Local documentation links must resolve inside this tree; proposed implementation files and future reports are not existing review dependencies.

Specifications define inputs/outputs, units, equations, precision, boundaries, ownership, errors and independent examples. Plugin specifications additionally define registration order, types, defaults, normalization and the supported feature subset. State each required behavior in its corresponding operator or interface specification.

Every phase closes only after its supported functions pass old/new comparisons through the public VS interface, including scalar and Highway paths. Core tests alone do not close a phase. AVS matching belongs to the final batch.
