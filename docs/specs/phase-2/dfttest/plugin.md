# DFTTest: phase-2 VS interface addendum

Specification: DFT-VS-002. Registration, argument order, types and defaults are unchanged from [DFT-VS-001](../../phase-1/dfttest/plugin.md). Apply that document except for these explicit overrides:

| Parameter | Default | Phase-2 rule |
| --- | --- | --- |
| tbsize | int 3 | Odd integer in 1..15 inclusive; also <=clip.num_frames. Valid set is 1,3,5,7,9,11,13,15 |
| tmode | int 0 | Only 0; reject 1 and other values at creation |
| tosize | int 0 | Parse int32, then force effective 0 for tmode=0 **before** overlap validation, even for supplied negative/large int32 values |
| twin | int 7 | Every ID 0..11, evaluated at actual length tbsize, including length 1 |
| tbeta | float 2.5 | Finite, >=0; used by temporal Kaiser, validated even if inactive |
| zmean | bool true | One window-shaped mean per full 3D block; [DFT-TIME-001](kernel-temporal.md) |

Defaults now enable T=3, so an otherwise valid clip with N>=3 succeeds without specifying tbsize. N=1 or 2 with omitted tbsize fails creation; explicitly tbsize=1 is supported. Never coerce an even T to odd, silently cap it at 15/N, reduce it at boundaries, or interpret it as a radius.

Both spatial modes and ftype=0..4 remain mandatory. Spatial smode=0 still requires odd sbsize; smode=1 may use even sbsize. Do not confuse spatial and temporal parity checks. Existing sosize normalization/divisibility, sigma2 default, numeric domains and exponent dispatch are unchanged.

## Deferred controls remain explicit

nlocation, slocation, ssx, ssy and sst must be empty; alpha must equal its phase-1 default (5 for ftype=0, otherwise 7); ssystem=0, dither=0 and existing execution restrictions remain. Even well-formed nonempty nlocation is an unsupported request, not a partially implemented profile. Non-finite values still fail when inactive. Noise estimation is phase 3 and contributes no phase-2 frame dependencies.

## Creation and evaluation

Parse/normalize and validate the full configuration before fetching frames, including with planes=[]. Check 3D volume/strides and all selected-plane geometry, build the window, stored-float energy/calibrated constants and optional G. A zero/non-finite mean-template DC with zmean=true is a creation error. Checked allocation failures are controlled failures.

For each n, request the [dependency set](../execution.md), retain its frames, gather T ordered slots, execute the [3D operator](kernel-temporal.md), and publish only n. Preserve properties from clip[n] and unchanged size/format/count/rate; copy unselected visible planes bitwise from n. Omitted planes selects all; explicit planes=[] selects none and needs only clip[n] after full argument validation. Unsupported parameters never become valid by selecting no planes.

Examples: tbsize=15,N=15 is legal; 0,2,4,14,16,17 and T>N fail at creation. tbsize=3,tosize=-123 is valid after normalization if other settings are valid. A supplied integer outside int32 still fails before normalization. Boundary replication is not a workaround for T>N.

Frame errors retain phase-1 owned-message/lifetime rules and must release workspace and frame references. No changed output or partial accumulator may escape a failed request.
