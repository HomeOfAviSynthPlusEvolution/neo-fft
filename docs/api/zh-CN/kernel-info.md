# KernelInfo

[目录](README.md)

查询当前进程可用的自动内核分派与 FFT 配置，不创建滤镜、不请求视频帧。

```python
import vapoursynth as vs
core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
info = core.neo_fft.KernelInfo()
print(info)
```

VS 调用无输入参数，返回 Python 字典，不返回 VideoNode，也不修改帧属性。

| 字段 | 类型 | 含义 |
|---|---|---|
| `target` | str | 自有频谱核自动分派目标，例如 AVX2、AVX3_SPR；具体字符串随构建和机器变化 |
| `fft_backend` | str | 实际 FFT 后端实现标签，例如 pocketfft-avx512；仅诊断，滤镜同名兼容输入的内容不影响它 |
| `fft` | str | FFT 配置标签，例如 pocketfft-native；注册为可选输出，当前实现总是提供 |
| `fft_lanes` | int | 所选 FFT 实现的 float SIMD lane 数，例如 16；注册为可选输出，当前实现总是提供 |

查询结果没有 `fft_threads` 或 `threads` 返回字段。DFTTest 的同名兼容输入参数不属于这个查询。`fft_lanes=16` 不表示 16 个线程。

查询报告自动路径，不接收实例或 opt。某实例使用 `opt=1` 时，自有核可以是标量，而查询的 target 仍是自动目标；FFT 又有独立分派。不要用它证明每一个块都运行某个专用 codelet，或据此推断速度。

需先显式加载插件或由宿主正常加载；不存在的插件命名空间、额外输入参数会报错。查询字符串是诊断标签，不应作为跨版本稳定枚举硬编码。字段选择原理见 [KernelInfo 原理](../../knowledge/zh-CN/kernel-info.md)。

## AviSynth

neo_fft_KernelInfo() 无参数，返回原生数组，固定顺序为 `[fft_backend, target, fft, fft_lanes]`；前三项为字符串，最后一项为整数。诊断含义和 VS 字段相同，没有线程计数字段。

```avs
LoadPlugin("/path/to/neo-fft.dll")
info=neo_fft_KernelInfo()
Assert(IsString(info[0]) && IsInt(info[3]))
return BlankClip(width=128,height=96,length=1,pixel_type="Y8")
```
