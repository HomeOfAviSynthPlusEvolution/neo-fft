# neo-fft 计算原理知识库

[目录](README.md)

这里沿像素到频谱、滤波再到输出的顺序解释当前实现。完整调用签名、参数默认值与错误处理见 [API 参考](../../api/zh-CN/README.md)。当前提供 VapourSynth 与 AviSynth 入口；两种滤镜内部均单线程，宿主仍可并发请求帧。

## 从滤镜开始

| 文章 | 内容 |
|---|---|
| [FFT3D](fft3d.md) | 空间/时间 Wiener、ROI、场排列、Kalman 与增强的处理链；4×4 完整例子 |
| [DFTTest](dfttest.md) | 时空体积、去均值、频谱增益、中心输出与重叠合成；完整算例 |
| [KernelInfo](kernel-info.md) | 自动分派、FFT 配置、lane 与线程的区别，以及诊断范围 |

## 公共计算基础

| 文章 | 解决的问题 |
|---|---|
| [样本域](shared/sample-domain.md) | 位深、float、色度中点、sigma 尺度和量化位置 |
| [块几何](shared/block-geometry.md) | 两算法不同的块网格、padding、反射和有效区域 |
| [FFT 约定](shared/fft-conventions.md) | 指数符号、半谱布局、DC/Nyquist、功率与归一化 |
| [窗口与重建](shared/windows-reconstruction.md) | FFT3D 成对肩窗、DFTTest 12 种窗、窗能量和 overlap-add |
| [执行与精度](shared/execution-precision.md) | 宿主并发、工作区、累加顺序、有限性与比较范围 |

## FFT3D 专题

| 文章 | 内容 |
|---|---|
| [噪声与增强](fft3d/noise-enhancement.md) | 均匀/曲线/采样优先级，degrid、Wiener、锐化、去光晕及预览 |
| [Kalman 状态与回放](fft3d/kalman.md) | 两步递推、初始化、重置、检查点、随机访问与状态内存 |
| [行频谱缓存](fft3d/spectra-cache.md) | 原始频谱、128 MiB 默认、登记淘汰、并发共享及 1080p/4K 算例 |

## DFTTest 专题

| 文章 | 内容 |
|---|---|
| [频谱模型](dfttest/spectral-models.md) | ftype=0..4、曲线与采样、sigma 标定、zmean 和阈值算例 |
| [时间重叠合成](dfttest/temporal-ola.md) | 绝对时间网格、偶数块、端点夹取、窗口归一化与累加顺序 |
| [Dither](dfttest/dither.md) | 整数输出量化、误差扩散、随机扰动、seed 和访问顺序确定性 |

使用者可先读滤镜主文，再按所用模式查专题。维护者建议先读样本域 → FFT 约定 → 窗口与重建 → 块几何 → 滤镜主文 → 执行与精度。所有公式例子都区分理想算术与实际 binary32；缓存 miss 比例不等同于总速度提升。

当前没有 DFTTest 跨帧频谱缓存，也没有可切换 FFT 后端。threads 等旧执行参数仅由宿主接受，插件完全忽略；兼容规则见 [API](../../api/zh-CN/README.md#线程与执行参数)。替换旧 Neo 系时需要修改的调用及结果差异见[迁移说明](../../api/zh-CN/migration.md)。
