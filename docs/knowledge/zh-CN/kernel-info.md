# 怎样理解内核与 FFT 分派

[目录](README.md)

KernelInfo 描述当前进程的自动选择，不是某个滤镜实例的运行追踪。调用和字段类型见 [API](../../api/zh-CN/kernel-info.md)。

## 版本与二进制身份

项目版本以 CMakeLists.txt 的 project VERSION 为唯一来源，构建时生成宿主共用头文件和 Windows DLL 版本资源。VS 数值注册版本只有主、次版本；完整三段版本可从插件描述或 DLL 产品版本读取。版本号并不区分同版本下的本地修改；复现问题时还应记录源码 commit、构建配置和二进制校验和。KernelInfo 的 SIMD/FFT 字段不能代替这些身份信息。

## 两套独立选择

自有核用 Highway 对 CPU 和构建能力作运行时分派；FFT 有独立的 PocketFFT 实现配置。FFT3D 的 opt=1、DFTTest 的 opt=1 只选择自有标量算子。其他受支持 opt 值使用自动路径，不代表数字对应的固定 ISA。

所以 target 与 fft_backend 的名字不必相同，也不能从 target 推断 FFT 的 lane 数。`fft` 表示配置选择，`fft_backend` 表示实际实现，`fft_lanes` 是该 FFT 的 float SIMD 宽度。

## 读一组结果

例如某次查询返回 `target=AVX3_SPR`、`fft=pocketfft-native`、`fft_backend=pocketfft-avx512`、`fft_lanes=16`。含义是自有自动核选择该 Highway 目标，native FFT 配置选中 AVX512 实现，每个向量有 16 个 float lane。它不是 16 个线程。

同一进程里创建 `opt=1` 的 DFTTest，并不会改变查询的自动 target；实例自有核使用标量，FFT 仍可 SIMD。CPU 不支持或构建未包含相应目标时使用可用回退，不执行不支持指令。

## 专用块优化的范围

DFTTest 的部分线性 Wiener 路径会按固定组批量 FFT，一些体积尺寸有专用变换；tmode=0 的合格路径可仅求需要的逆时间中心。其他类型、指数、尺寸或 OLA 使用对应回退。KernelInfo 不能证明某次调用进入哪一个 codelet，也不能把 lane 数直接当作加速倍数。

两个滤镜及 FFT 内部均单线程；threads、fft_threads、fft_backend 等兼容输入完全忽略，KernelInfo 不返回线程数字。同名 fft_backend 只读诊断报告实际实现，不受兼容输入影响。具体精度保证见[执行精度](shared/execution-precision.md)。
