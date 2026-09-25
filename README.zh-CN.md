# neo-fft

[English](README.md) | **简体中文** | [日本語](README.ja.md)

neo-fft 是 VapourSynth 和 AviSynth 的频域滤波插件，将 Neo FFT3D 与 Neo DFTTest 独立重写为一个插件，提供空间与时间降噪、Kalman 递推、频域锐化、去光晕及可配置的频谱滤波。

内部使用 C++17，提供标量内核、基于 Google Highway 的跨平台 SIMD，以及基于 PocketFFT 的变换实现。DualSynth2 将同一计算核心连接到两个宿主。VapourSynth 使用 `core.neo_fft`，AviSynth 使用带有 `neo_fft_` 前缀的函数。

## 设计

neo-fft 将频域计算与宿主帧管理分离。核心负责分块、加窗、FFT、频谱滤波及重建，宿主层负责参数、帧请求、属性和输出分配。FFT3D 与 DFTTest 共享基础设施，各自保留窗口、噪声模型和重叠合成规则；核心可以独立构建和测试。

实现依据行为规格开发，以标量实现和独立数学定义核对计算，再验证 SIMD 路径及固定参考插件的公开行为。相同参数不保证与所有历史 FFT3DFilter、DFTTest 或 Neo 版本逐像素一致。浮点舍入、阈值分支、Kalman 状态策略和 dither 均可能影响结果，具体差异见[迁移说明](docs/api/zh-CN/migration.md)。

两个滤镜和 FFT 内部均在宿主调用线程上计算，不创建工作线程或线程池。宿主仍可并发请求多个帧；SIMD 和批量 FFT 不表示内部多线程。

## 支持的操作

| 函数 | 用途 |
|---|---|
| `FFT3D` | 单帧或 2–5 帧 Wiener 降噪、Kalman 递推、锐化和去光晕，支持频率相关及采样噪声模型。 |
| `DFTTest` | 五种频谱滤波类型、空间与时间重叠合成、十二种窗口、频率曲线、噪声采样及 8 位输出 dither。 |
| `KernelInfo` | 查询自动选择的 SIMD 目标、FFT 实现及向量宽度。 |

两个视频滤镜支持固定格式、固定尺寸的平面 GRAY/YUV/RGB，样本类型为 8/10/12/14/16 位整数或 32 位浮点。AviSynth 另支持平面 YUVA/RGBA。输出保持输入的尺寸、格式、帧数和帧率；块尺寸及边界要求按每个实际处理的平面检查。

默认处理全部非 Alpha 平面，Alpha 默认复制。`planes=[0]` 只处理第一个平面，`planes=[]` 不处理任何平面。AviSynth 也接受旧式 `y/u/v/a` 平面模式；其中模式 1 明确表示不写入该输出平面，适用于随后丢弃该平面的脚本。具体优先级见 API。

FFT3D 默认 `bt=3`，使用相邻三帧；`bt=1` 为纯空间降噪，`bt=0` 为 Kalman。Kalman 默认最多预热 8 个历史帧，并可复用附近检查点，避免跳转时从头回放整个视频；缓存和请求历史可能影响递推结果。普通时间滤波的原始频谱缓存默认预算为 128 MiB，可通过 `cache_mb` 和 `cache_frames` 调整；这不是整个进程的内存上限。

DFTTest 默认 `tbsize=1`，仅使用当前帧；增加 `tbsize` 可启用时间滤波，`tmode` 控制中心输出或时间重叠合成。两个滤镜都不估计运动向量，也不做运动补偿。

## 文档与使用

API 文档说明怎样调用函数，知识库说明输入怎样经过窗口、变换、频谱模型和重建成为输出。

- [API 使用参考](docs/api/zh-CN/README.md)：函数签名、参数、默认值及使用示例。
- [计算原理知识库](docs/knowledge/zh-CN/README.md)：数据表示、公式、计算顺序、边界及精度。
- [从旧 Neo 系滤镜迁移](docs/api/zh-CN/migration.md)：必须修改的调用与结果不一致项。

可以显式加载构建出的插件，也可以将其放入 VapourSynth 的插件自动加载目录。下面使用 Windows 文件名；Linux 使用 `neo-fft.so`，其他平台请替换为实际插件路径。VapourSynth 插件标识符为 `org.neofilters.neo_fft`。

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")

print(core.neo_fft.KernelInfo())
clip = core.std.BlankClip(width=640, height=360, format=vs.YUV420P8, length=24)
output = core.neo_fft.FFT3D(clip, sigma=2.0, bt=3, planes=[0])
# 另一种选择：对原输入使用 DFTTest 空间降噪。
# output = core.neo_fft.DFTTest(clip, sigma=8.0, tbsize=1, planes=[0])
output.set_output()
```

这个最小示例使用合成剪辑展示调用方式。两个滤镜的 `sigma` 有不同的数学含义，不能通过填写相同数值得到相同降噪强度；实际使用时应分别调整。

同一个插件文件也提供 AviSynth C++ 接口，需要支持接口版本 11 的宿主，通过 `LoadPlugin` 加载：

```avs
LoadPlugin("/path/to/neo-fft.dll")
clip = BlankClip(width=640, height=360, length=24, pixel_type="YV12")
return neo_fft_FFT3D(clip, sigma=2.0, bt=3, planes=[0])
# DFTTest 的替代调用：
# return neo_fft_DFTTest(clip, sigma=8.0, tbsize=1, planes=[0])
```

AviSynth 数组使用 `[0, 1]` 等原生语法，数组参数也接受单值简写。DFTTest 的曲线及噪声采样位置另接受数字字符串，分隔符可以是空白、逗号或冒号。两个滤镜转发输入音频和奇偶性，输出帧属性来自对应的源帧。

显式表达参数省略语义时，VapourSynth 使用 `None`，AviSynth 使用 `Undefined()`。例如 `planes=None` 或 `planes=Undefined()` 使用默认平面选择，后者也允许 AVS 的 `y/u/v/a` 生效；空数组 `[]` 则明确选择零个平面。

FFT3D 的 `mt/ncpu/measure`、DFTTest 的 `threads/fft_threads` 及两者的 `fft_backend` 仅作为兼容输入接受，内容完全忽略；它们不会创建线程或切换 FFT 后端。旧脚本建议改用命名参数，参数顺序和迁移限制见对应 API。

## SIMD 与 CPU 选择

启用 SIMD 的构建会自动选择当前 CPU 支持且已编译的 Highway 目标，并保留标量回退。FFT 有独立的实现分派，部分固定尺寸还使用专用变换路径。

两个滤镜的 `opt=1` 选择自有标量内核，但不强制 FFT 为标量。其余受支持值走自动 SIMD，不能用旧 `opt` 数值指定或限制 ISA。需要整个构建关闭 SIMD 时，使用 `NEO_FFT_ENABLE_SIMD=OFF`。

`core.neo_fft.KernelInfo()` 返回包含 `fft_backend`、`target`、`fft` 和 `fft_lanes` 的字典，按字段名读取。AviSynth 的 `neo_fft_KernelInfo()` 返回数组，固定顺序为 `[fft_backend, target, fft, fft_lanes]`。查询描述自动路径，不是某个滤镜实例的执行追踪；`fft_lanes` 是向量通道数，不是线程数或加速倍数。

FFT 与自有内核的目标可能不同，更宽的 SIMD 也不保证更快。详见 [KernelInfo](docs/knowledge/zh-CN/kernel-info.md) 和[执行与精度](docs/knowledge/zh-CN/shared/execution-precision.md)。

## 构建与测试

需要 CMake 3.24 或更新版本、Git 及支持 C++17 的编译器。CMake 获取固定版本的 DualSynth2、PocketFFT，启用 SIMD 时还会获取 Highway 1.4.0。两个宿主的 SDK 均可从本地发现或自动下载。

下面构建双宿主插件和核心测试，不要求本地已安装视频宿主。Windows 默认包含 AVS C++ 入口，需要 MSVC 或 clang-cl；MinGW 构建须关闭 AVS 入口。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DNEO_FFT_TEST_VAPOURSYNTH=OFF
cmake --build build/release --config Release --parallel 4
ctest --test-dir build/release -C Release --output-on-failure
```

| 选项 | 用途 |
|---|---|
| `NEO_FFT_BUILD_VAPOURSYNTH=OFF` | 禁用 VapourSynth 入口。 |
| `NEO_FFT_BUILD_AVISYNTH=OFF` | 禁用 AviSynth 入口；两个宿主选项均为 OFF 时只构建核心。 |
| `NEO_FFT_ENABLE_SIMD=OFF` | 禁用 SIMD 内核和向量化 FFT 配置。 |
| `BUILD_TESTING=OFF` | 不构建测试。 |
| `NEO_FFT_TEST_VAPOURSYNTH=OFF` | 保留核心测试和插件，跳过 VapourSynth 宿主测试；启用 VS 入口和测试时默认 ON。 |
| `Python3_EXECUTABLE=/path/to/python` | 指定宿主测试所用 Python；VS 测试需该环境可导入 VapourSynth 并加载匹配架构的运行时。 |
| `NEO_FFT_VS_SDK=/path/to/sdk` | 指定本地 VapourSynth SDK。 |
| `NEO_FFT_AVS_SDK=/path/to/sdk` | 指定本地 AviSynth SDK。 |
| `NEO_FFT_TEST_AVISYNTH=ON` | 启用 AviSynth 宿主测试，默认关闭；需要 Python 和匹配架构的运行时。 |
| `NEO_FFT_AVISYNTH_RUNTIME=/path/to/avisynth.dll` | 指定 AviSynth 宿主测试使用的运行时库。 |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | 使用本地 DualSynth2 源码代替固定版本下载。 |

插件构建目标为 `neo_fft`，输出文件基本名为 `neo-fft`，默认包含两个宿主入口。测试覆盖 FFT、窗口、滤波算子、边界、标量/SIMD 比较及宿主行为。与旧插件比较的黑盒测试另需固定参考二进制。

CI 包含 Windows x64、Linux x64、macOS ARM64 及 Linux ASan/UBSan 检查。Linux runner 使用 Ubuntu 26.04，GCC 使用系统默认版本，Sanitizer 检查使用 Clang 22。发布工作流构建 Windows、Linux 和 macOS 的 x64/ARM64 产物，VapourSynth 宿主测试目前在 Windows x64 上运行，AviSynth 宿主测试通过上述选项单独启用。工作流产物的实际构建环境和测试范围记录在包内，不代表适用于所有 Linux 发行版。

## 性能

已有测量中，FFT3D 常用路径的吞吐量约为旧 Neo FFT3D 的 **1.72–3.05 倍**，DFTTest 常用 Wiener 路径约为旧 Neo DFTTest 的 **4.58–7.04 倍**。比值为 **neo-fft 吞吐量 / 参考滤镜吞吐量**，等价于参考滤镜耗时 / neo-fft 耗时；**大于 1 表示 neo-fft 更快**。下列范围跨越 8/16 位整数、32 位浮点及 AVX2/AVX-512 配置，不是置信区间。

| 常用路径 | 相对吞吐量 |
|---|---:|
| FFT3D 空间降噪（`bt=1`） | 2.39–2.85× |
| FFT3D 两帧降噪（`bt=2`） | 2.32–3.05× |
| FFT3D 三帧降噪（`bt=3`，默认时间模式） | 1.97–2.53× |
| FFT3D 四帧降噪（`bt=4`） | 1.86–2.37× |
| FFT3D 五帧降噪（`bt=5`） | 1.72–2.15× |
| FFT3D Kalman（`bt=0`，顺序请求） | 1.90–2.43× |
| DFTTest 空间 Wiener（`tbsize=1`） | 4.60–6.65× |
| DFTTest 三帧 Wiener（`tbsize=3, tmode=0`） | 5.25–7.04× |
| DFTTest 五帧 Wiener（`tbsize=5, tmode=0`） | 4.58–5.77× |

这些范围来自 AviSynth 下的单线程历史对照，参考为旧 Neo FFT3D / Neo DFTTest 的 Highway 现代化版本。测量不包含后续优化，也不代表整条处理链的吞吐量；实际结果随输入、参数、硬件和宿主并发数变化。

## 开发与贡献

维护者负责技术方向、变更审核和发布。欢迎问题报告、建议与贡献；修改数值语义、公开接口或重要架构前，建议先讨论目标和方案。

本项目使用 AI 辅助实现、测试和审查。贡献应说明问题、方案、验证方法和 AI 参与方式。报告问题请提供版本、系统、CPU、编译器、构建选项、输入输出格式及最小复现；数值差异还应注明参考版本、参数和请求顺序，性能报告应注明计时范围及线程配置。

## 致谢与许可证

感谢以下上游项目的作者与贡献者，他们的工作为 neo-fft 的接口与频域滤波功能提供了基础：

- [FFT3DFilter](https://github.com/pinterf/fft3dfilter)：由 Alexander G. Balakhnin（Fizick）最初开发，martin53 完成 AviSynth 2.6 适配，Ferenc Pintér（pinterf）继续改进并加入高位深支持。
- [DFTTest](https://github.com/pinterf/dfttest)：由 tritical 最初开发，Firesledge 加入 16 位处理，DJATOM 完成 AviSynth+ 移植，pinterf 继续完善高位深支持和跨平台构建。

neo-fft 还使用了以下基础库：

- [Google Highway](https://github.com/google/highway)：提供跨平台 SIMD 支持。
- [PocketFFT](https://github.com/mreineck/pocketfft)：用于 FFT3D 与 DFTTest 的频域变换。
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2)：连接 VapourSynth、AviSynth 与共享计算核心。

感谢参与测试、报告问题和改进的开发者与用户。

感谢 [烧饼论坛](https://sb.sb) 赞助本项目开发使用的 LLM 订阅。

neo-fft 采用 GNU 通用公共许可证第 2 版或更新版本（`GPL-2.0-or-later`），完整条款见 [LICENSE](LICENSE)。第三方组件保留各自的版权声明和许可条款。
