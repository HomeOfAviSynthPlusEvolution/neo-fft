# DFTTest phase-4 interface

Specification: DFT-VS-004. Extends [phase-3 interface](../../phase-3/dfttest/plugin.md). Registration order/defaults are unchanged except the now-active omitted dither seed has deterministic semantics.

| Parameter | Domain/default | Effect |
| --- | --- | --- |
| dither | int32, default 0; >=0 | 0 ordinary conversion, 1 diffusion, >=2 diffusion plus coordinate-based noise |
| dither_seed | Optional nonnegative int32 | Omitted resolves to 0; used only for dither>=2 on selected UInt8 planes |
| opt | int32, default 0 | Accept 0,1,2,3,8; execution mapping below |
| threads,fft_threads | Existing int32/default 0 | Positive requests now supported under [execution limits](../execution.md) |

Only selected **UInt8** output planes use [the dither operator](kernel-dither.md). UInt10/12/14/16 and Float32 keep ordinary conversion for all valid dither values; do not introduce high-bit-depth dithering. Mode 0 or 1 never initializes/advances a hidden RNG; seed has no effect. Negative dither/seed fails creation, including when inactive or planes=[]. There is no arbitrary small maximum dither mode; perform the specified binary32 scale conversion without integer overflow.

Dither runs after the complete temporal/spatial filter, inverse compensation, overlap-add/center extraction and crop. It does not modify the frequency-domain PSD, mean removal/restoration, filtering, windows, sample estimation or reconstruction. An unselected plane is a bitwise copy. planes=[] validates configuration then copies n without model, FFT or dither allocations.

Pixel coordinates for noise and diffusion are the **visible output plane's** coordinates, with x=y=0 at that plane's origin. n is the output frame index and p is the actual source/output plane index, not its position in the planes argument. No dependence on block, strip, worker, request order or retry count is allowed. Existing temporal and sampled-noise dependency closures remain unchanged.

tmode=1 remains unsupported, even with planes=[]. Dither does not open a new temporal output mode. Error metadata/ownership, source properties and output format/length remain inherited.
