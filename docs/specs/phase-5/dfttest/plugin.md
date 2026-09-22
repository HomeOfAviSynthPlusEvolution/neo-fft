# DFTTest phase-5 interface

Specification: DFT-VS-005. Extends [phase 4](../../phase-4/dfttest/plugin.md). Signature, argument order and defaults are unchanged. This document supersedes earlier tmode=1 rejection and unconditional odd-T restrictions only where stated.

| Control | Effective domain |
| --- | --- |
| tmode | int32, default 0; accept 0 or 1 |
| tbsize=T | int32, default 3; 1<=T<=15 and T<=N; odd for tmode=0, either parity for tmode=1 |
| tosize=O | int32, default 0; tmode=0 normalizes to 0 before overlap validation; tmode=1 requires 0<=O<T |
| Heavy temporal overlap | For tmode=1, if O>floor(T/2), require T divisible by T-O |

tmode=0 keeps its existing centered-window slot mapping, raw temporal window and center-only output. tmode=1 uses the [block lattice](kernel-temporal-geometry.md) and [normalized temporal window](kernel-temporal-ola.md). A time block produces a contribution at each of its T logical positions; output n combines every covering block's corresponding slice. O=0 is a valid non-overlapping block mode; do not turn it into centered mode. Omitted tosize remains 0, even for tmode=1.

Examples: T=4,O=2 and T=5,O=2 are valid; T=6,O=4 is valid heavy overlap; T=5,O=3 fails because 5 is not divisible by 2. T=4,tmode=0 fails. T=1,tmode=1 requires O=0 and still evaluates/normalizes its temporal window. It need not be bit-identical to tmode=0,T=1 for all window choices.

Validate types, integer ranges, raw finite values, mode/overlap rules and inherited array/sample-rectangle structure before derived arithmetic. planes=[] still performs this validation but returns a bitwise copy of source n with its properties and no FFT/window/model/cache work. Omitted planes selects all, explicit empty selects none. No invalid configuration becomes valid by selecting no planes.

All twelve temporal windows remain available subject to finite nonzero phase energy; zmean=true also requires a finite nonzero template DC. sbsize/smode/sosize constraints and per-plane spatial reflection admission remain unchanged. Temporal endpoint duplication is allowed independently of the spatial one-reflection restriction.

ftype=0..4, zmean, scalar/curve/nlocation precedence, sample amplitudes, signed float chroma, finite-content scope and output quantization remain inherited. Active nlocation tuples use exactly T chronological frames fn..fn+T-1, even when T is even; their first frame need not lie on the output block lattice. Validate fn+T<=N with checked arithmetic and do not clamp sample intervals. Curves support even T with the inherited full-axis frequency mapping.

Apply dither exactly once to the fully combined visible output; its seed/frame/plane/coordinate identity is unchanged. Unselected planes and properties always come from source n. Size, format, length and rate are unchanged; no multi-output clip or field-rate conversion is introduced.

Creation errors include invalid mode/parity/overlap, unsafe geometry, non-finite active tables and checked-size overflow. Source/model/transform overflow during processing fails that frame atomically. Do not fall back to tmode=0 or reduce O/T. Execution parameters retain phase-4 mappings.
