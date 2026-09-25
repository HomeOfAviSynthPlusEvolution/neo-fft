# 从旧 Neo 系迁移

[目录](README.md)

本文面向旧 **Neo FFT3D**（neo-fft3d）和 **Neo DFTTest**（neo-dfttest）的 AviSynth / VapourSynth 脚本。两者在 neo-fft 中合并为同一插件，保留独立算法入口。本文不承诺传统 FFT3DFilter、tritical/pinterf DFTTest 的全部历史接口兼容，例如 nfile、nstring、sstring 不属于这里的同名参数迁移范围。

常用 CPU 功能可以迁移，但“调用成功”和“与旧实现逐像素一致”是两个不同条件。以下先列必须修改项，再列可能改变结果的行为。完整参数见 [FFT3D](fft3d.md)、[DFTTest](dfttest.md) 和 [KernelInfo](kernel-info.md)。

## 必须修改项

只修改脚本实际涉及的项目；显式加载路径、函数名或命名空间需要换成新入口。

| 场景 | 旧调用 | 新调用或必要操作 |
|---|---|---|
| 显式加载插件 | neo-fft3d.dll、neo-dfttest.dll | 加载 neo-fft.dll；Linux 产物为 neo-fft.so。自动加载时安装新插件。 |
| AVS FFT3D | `neo_fft3d(c, ...)` | `neo_fft_FFT3D(c, ...)` |
| AVS DFTTest | `neo_dfttest(c, ...)` | `neo_fft_DFTTest(c, ...)` |
| VS FFT3D | `core.neo_fft3d.FFT3D(c, ...)` | `core.neo_fft.FFT3D(c, ...)` |
| VS DFTTest | `core.neo_dfttest.DFTTest(c, ...)` | `core.neo_fft.DFTTest(c, ...)` |
| 按插件标识符查找的工具 | 旧插件标识符 | 使用 `org.neofilters.neo_fft`。 |
| clip 之后使用位置参数 | 依赖旧签名顺序 | 按新 API 顺序改写，建议改成命名参数。新参数表不保证与旧参数表逐槽对应。 |
| AVS DFTTest 的 YUVA/RGBA 输入，希望保留旧默认 Alpha 处理 | 省略 `a`，旧默认为 3 | 没有显式 `planes` 时添加 `a=3`；已有 `planes` 时将索引 3 加入其中。 |
| 依赖旧诊断字段或返回结构 | 读取线程数等旧字段 | 按新 KernelInfo 接口改写；VS 返回字典，AVS 返回四元素数组，不再返回线程数字段。 |

AVS 需要支持接口版本 11 的 AviSynth+ 宿主；较旧宿主需要升级。Windows AVS 构建使用 MSVC 或 clang-cl。两种宿主入口的构建与加载说明见 [API 目录](README.md)。

以下以已有输入 `c` 为例，分别将旧函数名换成新函数名，参数使用名称传入：

```avs
LoadPlugin("/path/to/neo-fft.dll")
# 旧调用：neo_fft3d(c, sigma=2.0, bt=3, y=3, u=2, v=2)
return neo_fft_FFT3D(c, sigma=2.0, bt=3, y=3, u=2, v=2)
```

```python
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
# 旧调用：core.neo_dfttest.DFTTest(c, sigma=8.0, tbsize=3, planes=[0])
output = core.neo_fft.DFTTest(c, sigma=8.0, tbsize=3, planes=[0])
```

### 可以保留的参数，以及不再生效的控制

AVS 的 `y/u/v/a` 和 VS 的 `planes` 均可使用；AVS 也接受 `planes`。显式 `planes` 优先于全部 `y/u/v/a`，包括显式空数组；两个滤镜的 `planes=[]` 均全部复制。VS `planes=None` 和 AVS `planes=Undefined()`（必须带括号）等价于省略，不覆盖 AVS 的 y/u/v/a。未指定选择参数时，两个滤镜均处理现有非 Alpha 平面并复制 Alpha。没有 `planes` 时，未提供的模式分别默认 `y=u=v=3, a=2`；1 跳过写入，2 复制，3 处理。模式 1 与旧接口一致，不保证该输出平面的像素内容；可继续用于只取有效平面再组合的脚本。RGB 中对应 R/G/B/A。

已验证的旧 Neo FFT3D r14 和 Neo DFTTest r12 的 VS 接口均拒绝显式 `planes=[]`；当前允许空选择属于扩展，不是这些旧脚本必须修改的项目。先前 neo-fft 开发版本曾将 FFT3D 的 `[]` 当作默认选择；依赖该行为的脚本请改为省略参数或 VS `None` / AVS `Undefined()`。

DFTTest 的 `nlocation`、`slocation`、`ssx`、`ssy`、`sst` 在 AVS 下接受空白、逗号、冒号及混合分隔的字符串，也接受原生数字数组；不需要为了分隔符重写旧字符串。VS 继续传数字数组。小数点固定为 `.`，逗号始终分隔数值。例如 `slocation="0:4, 1:12"` 等价于 `slocation=[0,4,1,12]`。

FFT3D 的 `mt/ncpu/measure`、DFTTest 的 `threads/fft_threads` 和两者的 `fft_backend` 可以保留，宿主检查注册类型，但插件不读取其内容。它们不创建内部线程，也不切换 FFT 后端，内部固定使用 PocketFFT。若旧脚本依赖这些参数分配并行度，需要改为配置宿主的帧级并发。

`opt=1` 选择本项目标量核，其他允许值使用自动 SIMD；这些值不再指定或限制某个 ISA，`opt=1` 也不强制 FFT 实现为标量。FFT3D 接受任意 int32，DFTTest 接受 0/1/2/3/8。依赖旧 opt 数值限制 ISA 的脚本不能继续把它当作相同控制；详见各滤镜的 API。

## 结果不一致项

以下既包含有意定义的行为，也包含已有参考对照中的数值例外；不能统一归为“只有 FFT 舍入误差”。参考差异针对已验证的旧构建，不代表所有历史发布件都有相同行为。

| 场景 | 当前行为及迁移影响 |
|---|---|
| AVS Alpha 默认选择 | 默认复制 Alpha；旧 Neo DFTTest 默认处理。要保留处理行为须显式选择，见上表。Alpha 按全分辨率、亮度/RGB 样本域处理。 |
| FFT3D 普通 Wiener 的正 `pfactor` | 按强度线性缩放采样噪声功率；固定旧参考只区分零和非零。非 1 的正值可能改变降噪强度。对照时用双方 `pfactor=1`；这不是让用户放弃其他强度。Kalman 是例外，正值只选择采样，不缩放功率。 |
| FFT3D `pshow` | 使用独立预览窗口和实际平面的色度标记，避免旧参考的窗口状态和亮度偏移问题。重复、乱序预览及部分平面的结果可能不同；不复刻这些旧行为。 |
| FFT3D Kalman | 使用 kalman_warmup=8 的有界预热与附近已完成检查点。冷跳转不再回放全部前缀，缓存、请求历史和并发可能改变结果；既不同于旧实现的调用历史策略，也不同于先前 neo-fft 的规范全历史回放。即使正序请求，阈值附近舍入也可能改变重置分支，详见下文。 |
| DFTTest `dither>=2` | 随机扰动由 seed、帧号、平面和像素坐标确定，不追随旧实现随工作槽或请求顺序演进的随机序列。同一个 seed 不保证产生旧噪声图案。 |
| 浮点计算与量化 | FFT 实现、SIMD 和运算分组可能改变舍入；整数也可能因阈值或最终量化而不同。历史测试容差只覆盖相应输入和配置，不是任意输入的统一误差上限。 |
| 非法或退化配置 | 非有限值、重复曲线节点、越界采样等会被明确拒绝；算法参数即使在当前模式下不参与计算，也可能仍须通过原始校验。不会复刻越界或除零行为；显式平面模式 1 的不写入语义保留。完全忽略的执行兼容参数不适用这些算法值校验。 |

### Kalman 已知参考例外

已有固定参考 `bc9dd394` 的 AVS 对照中，6 个配置（B32、sigma=2、kratio=2、128×96 RGB16/float，覆盖三种派发）超出了预先确定的像素容差：RGB16 最大差异 33，原预算 4；float 最大约 7.9751e-4，原预算 2e-5。它们仍记录为参考不一致，没有放宽容差，也不能计作通过。

已定位到加窗/FFT 舍入使严格运动阈值分支分歧，随后递推和重建放大差异；未发现当前实现违背既定算子规则。这既不是“所有 Kalman 输入都有这个误差”，也不是“已证明所有旧 Kalman 实现错误”。需要保持既有逐像素金样的使用者应单独检查所用 Kalman 配置。

计算背景见同语言知识库：[噪声与增强](../../knowledge/zh-CN/fft3d/noise-enhancement.md)、[Kalman 状态与回放](../../knowledge/zh-CN/fft3d/kalman.md)、[Dither](../../knowledge/zh-CN/dfttest/dither.md)及[执行与精度](../../knowledge/zh-CN/shared/execution-precision.md)。
