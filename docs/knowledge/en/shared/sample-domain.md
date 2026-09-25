# Sample units and output conversion

[Contents](../README.md)

Both filters use float32 work buffers, but their pixel values have different meanings. Here bits is integer input depth and q is amplitude before the analysis window. Neither filter automatically remaps video range to full range.

| Input | FFT3D q | DFTTest q |
|---|---|---|
| Integer luma, RGB or Alpha | sample | sample / 2^(bits−8) |
| Integer YUV chroma | sample − 2^(bits−1) | sample / 2^(bits−8), no midpoint subtraction |
| float32, all planes | sample | sample × 255 |

Integer luma values 128, 512 and 32768 at 8/10/16 bits all become 128 in DFTTest. FFT3D preserves their native scale and multiplies sigma by 1, 4 and 256. Float luma 0.5 becomes 127.5 in DFTTest but stays 0.5 in FFT3D, whose sigma is divided by 255. Neutral 10-bit chroma 512 becomes 0 in FFT3D and 128 in DFTTest. Float chroma −0.25 becomes −63.75 in DFTTest; neither subtracts another midpoint from float chroma.

FFT3D sigma is a standard deviation, squared to form noise power. DFTTest sigma is a type-dependent power or gain: do not square it or interpret equal numeric settings as equivalent between filters.

## Final clipping

Let z be the reconstructed value. FFT3D integer output evaluates `(z+0.5)+base`, safely truncates and clips to `[0,2^bits−1]`; only integer chroma has a nonzero base. Float output clips to `[0,1]`, including float chroma, so processed negative chroma may become 0.

DFTTest ordinary UInt8 conversion is `clip(trunc(z+0.5),0,255)`; higher integer depths use `clip(trunc(z*2^(bits−8)+0.5),0,2^bits−1)`; float uses `z/255` without clipping. Optional UInt8 [dither](../dfttest/dither.md) replaces ordinary conversion.

For example, DFTTest z=128.25 produces 128 at u8, 513 at 10 bits and approximately 0.502941 as float. Quantizing to u8 before restoring bit depth would lose precision. Saturation precedes unsafe integer conversion; intermediate blocks are not clipped individually.

## Copies, formats and finite values

Both functions support fixed planar GRAY/YUV/RGB with one or three planes, plus four-plane YUVA/RGBA on AviSynth: integer 8/10/12/14/16 or float32. Duplicate planes are selected once. Both default to all existing non-Alpha planes. Explicit `planes=[]` selects none in both filters, copying everything after configuration validation. VS `planes=None` and AVS `planes=Undefined()` mean omission, not an empty list; AVS omission still uses y/u/v/a and their defaults. Explicit index 3 selects Alpha on AVS; its dimensions and sample scaling follow luma/RGB, never subsampled chroma. AVS y/u/v/a modes apply only without explicit planes, default to 3/3/3/2, and select skip-write for 1, copy for 2, or process for 3. Mode 1 builds no processing plan for that plane and skips writes in ordinary and Kalman output, including frame 0; output pixels are unspecified and the full output frame is still allocated. Other unselected planes, areas outside selected FFT3D ROIs and Kalman output frame 0 are bitwise copies.

Consumed float samples and active intermediates must be finite. DFTTest sampling may consume a patch from an output-unselected plane, which must also be finite. Copy regions do not fail merely for containing NaN. Dimensions, frame count and rate are unchanged; properties come from source n.

See [FFT3D API](../../../api/en/fft3d.md), [DFTTest API](../../../api/en/dfttest.md) and [execution precision](execution-precision.md).
