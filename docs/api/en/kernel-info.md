# KernelInfo

[Contents](README.md)

Reports automatic kernel dispatch and FFT configuration in the current process, without creating a filter or requesting frames.

```python
import vapoursynth as vs
core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
info = core.neo_fft.KernelInfo()
print(info)
```

The VS call takes no arguments and returns a Python dictionary, not a VideoNode. It does not change frame properties.

| Field | Type | Meaning |
|---|---|---|
| `target` | str | Automatic target for own spectral kernels, such as AVX2 or AVX3_SPR; actual strings depend on build/machine |
| `fft_backend` | str | Actual FFT implementation label, such as pocketfft-avx512; diagnostic only, unaffected by the ignored same-name compatibility input |
| `fft` | str | FFT profile label, such as pocketfft-native; registered optional, always supplied by the current implementation |
| `fft_lanes` | int | Float SIMD lane count of the selected FFT implementation, such as 16; registered optional, always supplied currently |

There is no `fft_threads` or `threads` output. DFTTest's same-name compatibility inputs are separate. `fft_lanes=16` does not mean 16 threads.

This query describes automatic dispatch, not a particular instance or opt setting. An `opt=1` instance can use scalar own kernels while target still reports the automatic target. FFT dispatch is separate again. This cannot prove that every block uses a specific codelet or predict speed.

Load the plugin explicitly or through normal host loading. A missing plugin namespace or extra arguments cause errors. Strings are diagnostic labels, not stable cross-version enums. See [dispatch interpretation](../../knowledge/en/kernel-info.md).

## AviSynth

neo_fft_KernelInfo() takes no arguments and returns a native array ordered `[fft_backend, target, fft, fft_lanes]`: three strings followed by an integer. Meanings match the VS fields, with no thread-count field.

```avs
LoadPlugin("/path/to/neo-fft.dll")
info=neo_fft_KernelInfo()
Assert(IsString(info[0]) && IsInt(info[3]))
return BlankClip(width=128,height=96,length=1,pixel_type="Y8")
```
