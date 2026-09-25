# neo-fft API 参考

[目录](README.md)

这里说明公开函数怎样调用：输入格式、参数类型与默认值、返回结果、可运行示例和常见错误。计算公式、处理顺序与数值例子见[计算原理知识库](../../knowledge/zh-CN/README.md)。

VapourSynth 使用命名空间 `core.neo_fft`，插件标识符为 `org.neofilters.neo_fft`。AviSynth 使用 `neo_fft_` 函数前缀；两种宿主均提供 FFT3D、DFTTest 和 KernelInfo。

下列三篇覆盖当前全部公开入口；每篇包含完整参数、返回值、示例及限制。

| 函数 | 用途 |
|---|---|
| [FFT3D](fft3d.md) | 空间/时间频域降噪、Kalman 递推以及锐化和去光晕；返回 VideoNode |
| [DFTTest](dfttest.md) | 五类频谱滤波、时空重叠、曲线/采样与 dither；返回 VideoNode |
| [KernelInfo](kernel-info.md) | 当前自动 SIMD 与 FFT 配置；返回诊断字典 |

从旧 Neo FFT3D / Neo DFTTest 切换时，先读[迁移说明](migration.md)，区分必须修改的调用与结果差异。

## 版本标识

首版版本为 `0.9.0`。VS 的插件描述和 AVS 的加载描述包含完整版本；VS 的数值注册版本仅表示主、次版本（当前为 0.9），不包含补丁号。Windows DLL 的产品版本为 `0.9.0`，文件版本为 `0.9.0.0`。KernelInfo 保持原有字段，报告计算配置。

## 线程与执行参数

FFT3D 与 DFTTest 均在当前宿主调用线程内完成一个请求，不创建内部工作线程或线程池。PocketFFT 内部也固定为单线程。SIMD 与批量 FFT 仍然启用；它们不要求额外线程。宿主可并发请求多个输出帧，帧级并发由宿主控制。

| 接口 | 当前行为 | 脚本迁移 |
|---|---|---|
| DFTTest `threads`、`fft_threads` | 可选整数兼容参数，接受后完全忽略，不读取、归一化或保存 | 可保留原调用参数，不会创建线程 |
| FFT3D `mt`、`ncpu`、`measure` | 可选兼容参数，类型依次为布尔、整数、布尔；接受后完全忽略 | 可保留原调用参数 |
| FFT3D / DFTTest `fft_backend` | 可选字符串兼容参数，内容完全忽略；内部固定使用 PocketFFT | 可保留原调用参数，任何后端名称均不改变实际实现 |
| KernelInfo 返回值 `fft_backend` | 保留，报告实际 FFT 实现及 SIMD 路径 | 只读诊断，不用于选择后端 |
| KernelInfo 返回值 `fft_threads` | 已删除；保留 FFT 后端、配置与 SIMD 诊断字段 | 删除对该返回字段的读取；兼容输入参数不会新增返回字段 |

兼容参数仅登记在宿主签名中，不进入 DS2 参数解析或算法配置。宿主仍按上述类型检查传参，但插件不检查这些参数的数值范围或字符串内容；VS 整数兼容参数不受插件 int32 限制。新增兼容参数排在当前签名末尾，不移动已有参数；迁移旧位置参数调用时仍须核对顺序，建议按名称传入。

每个实例跨平面合计最多保留一个空闲工作区，且不超过 64 MiB；活动请求可独立取得工作区。这个空闲保留上限不限制宿主并发，也不是整个进程的内存上限。FFT3D 的原始频谱缓存另由 `cache_frames` / `cache_mb` 控制。

资源保留与请求调度见[执行与精度](../../knowledge/zh-CN/shared/execution-precision.md)。

## 加载插件

```python
import vapoursynth as vs
core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
print(core.neo_fft.KernelInfo())
```

Windows 产物为 neo-fft.dll，Linux 为 neo-fft.so；替换路径后即可使用。自动加载环境可省略 LoadPlugin。原理阅读从 [FFT3D](../../knowledge/zh-CN/fft3d.md) 或 [DFTTest](../../knowledge/zh-CN/dfttest.md) 开始。

## AviSynth 调用与构建

需要支持接口版本 11 的 AviSynth+ 宿主。使用 LoadPlugin 加载同一个 neo-fft.dll，然后调用 neo_fft_FFT3D、neo_fft_DFTTest、neo_fft_KernelInfo。共享参数的数值语义与 VS 一致。AVS 的 y/u/v/a 平面模式在 FFT3D 中位于 fft_backend 与 kalman_warmup 之间，在 DFTTest 中位于共享参数末尾；支持平面 YUVA/RGBA，Alpha 默认复制、可显式选择处理；数组写为原生 `[0,1]`，也接受单值作为一个元素的数组。两个滤镜的 `planes=[]` 均表示不处理任何平面，全部复制。省略参数默认处理全部现有非 Alpha 平面；AVS 仍可由 y/u/v/a 指定选择。显式表达省略语义时，VS 写 `planes=None`，AVS 写 `planes=Undefined()`（必须带括号）。省略值不覆盖 AVS y/u/v/a，而 `[]` 会覆盖。

两个视频滤镜转发输入音频和奇偶性，输出帧属性仍来自源帧 n。KernelInfo 返回原生数组，顺序见其 API 页。内部单线程；Prefetch 由宿主调度并发。

DFTTest 的采样位置与四个曲线参数还支持空白、逗号或冒号分隔的数字字符串，详见 [DFTTest 的 AviSynth 调用](dfttest.md#avisynth)。

构建默认同时启用 NEO_FFT_BUILD_AVISYNTH 和 NEO_FFT_BUILD_VAPOURSYNTH，可分别关闭以仅构建一个宿主；都关闭时只构建核心。Windows AVS C++ ABI 要求 MSVC 或 clang-cl，MinGW 构建需关闭 AVS。NEO_FFT_AVS_SDK 可指定 SDK；运行测试还需配置 NEO_FFT_TEST_AVISYNTH=ON 和 NEO_FFT_AVISYNTH_RUNTIME。产物 Windows 为 neo-fft.dll，Linux 为 neo-fft.so。
