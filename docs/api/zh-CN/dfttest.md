# DFTTest

[目录](README.md)

对同一位置的时空块做频域滤波，支持五种频谱增益、曲线和采样噪声模型，以及时间重叠合成。它不估计运动向量。

## 调用方式

```text
core.neo_fft.DFTTest(
    clip [, ftype, sigma, sigma2, pmin, pmax, sbsize, smode, sosize,
    tbsize, tmode, tosize, swin, twin, sbeta, tbeta, zmean, f0beta,
    nlocation, alpha, slocation, ssx, ssy, sst, ssystem, dither,
    dither_seed, planes, opt, threads, fft_threads, fft_backend]
) -> VideoNode
```

方括号表示可选参数；只有 clip 必填。参数顺序对应注册顺序，建议使用关键字。输入为固定格式、固定尺寸且至少一帧的平面 GRAY/YUV/RGB，1 或 3 平面（AviSynth 另支持四平面 YUVA/RGBA）；支持整数 8/10/12/14/16 位和 float32。AviSynth 对应函数为 neo_fft_DFTTest，共享参数顺序与此一致，末尾另加 AVS 专用的 y/u/v/a。

## 滤波参数

| 参数 | 类型 | 默认 | 取值和作用 |
|---|---|---|---|
| `clip` | VideoNode | 必填 | 输入剪辑 |
| `ftype` | int | 0 | 0 Wiener；1 硬阈值；2 直接增益；3 功率区间内外两种增益；4 平滑功率相关增益 |
| `sigma` | float | 8 | 非负；0/1 为功率参数，2/3/4 为增益。**不平方**，不能照搬 FFT3D 的标准差解释 |
| `sigma2` | float | 8 | 非负；类型 3 区间外增益，独立默认 8，不继承 sigma |
| `pmin`, `pmax` | float | 0, 500 | `0≤pmin≤pmax`；类型 3/4 使用，经窗能量标定 |
| `f0beta` | float | 1 | 大于 0，类型 0 增益的指数；其他类型不使用 |
| `zmean` | bool | True | 在整个时空体积上移除并恢复 DC 推导的窗形频谱；不是逐帧去均值 |

## 块与窗口

| 参数 | 类型 | 默认 | 取值和作用 |
|---|---|---|---|
| `sbsize` | int | 16 | 正整数，以每平面样本计；smode=0 时必须为奇数 |
| `smode` | int | 1 | 0：每个像素的中心块仅输出中心样本；1：空间重叠合成 |
| `sosize` | int | 12 | smode=0 时任何 int32 归一化为 0；smode=1 时 `0≤O<S`，若 `O>floor(S/2)`，要求 `S%(S-O)=0`，S=sbsize |
| `tbsize` | int | 3 | `1..15` 且不超过输入帧数；tmode=0 要求奇数，tmode=1 允许偶数 |
| `tmode` | int | 0 | 0：以当前帧为中心，只取时间中心输出；1：固定时间块网格，合成所有覆盖当前帧的块 |
| `tosize` | int | 0 | tmode=0 时任何 int32 归一化为 0；tmode=1 时 `0≤O<T`，若 `O>floor(T/2)`，要求 `T%(T-O)=0` |
| `swin`, `twin` | int | 0, 7 | 空间／时间原始窗，0..11；长度为 1 时也实际求值 |
| `sbeta`, `tbeta` | float | 2.5, 2.5 | 非负，对应 Kaiser 窗参数；其他窗不使用 |

窗编号：0 Hann、1 Hamming、2 Blackman、3 四项 Blackman–Harris、4 Kaiser、5 七项 Blackman–Harris、6 flat top、7 矩形、8 Bartlett、9 非对称 Bartlett–Hann 变体、10 Nuttall、11 Blackman–Nuttall。名称不能代替精确定义，见[窗口](../../knowledge/zh-CN/shared/windows-reconstruction.md)。所有活动窗口要求非零有限能量；zmean 另要求模板 DC 非零。

空间 padding 使用一次不重复端点的反射，各选中平面须满足反射范围；小色度平面也单独检查。时间越界使用首尾帧重复，保留完整 T 个逻辑位置，仍要求 T≤帧数。`tmode=1,T=4,O=2` 合法，`T=5,O=3` 不合法；省略 tosize 的默认值始终是 0。

## 曲线与采样

| 参数 | 类型 | 默认 | 取值和作用 |
|---|---|---|---|
| `slocation` | float[] | 空 | 共享频率曲线；非空时覆盖轴曲线 |
| `ssx`, `ssy`, `sst` | float[] | 各空 | X、Y、时间轴曲线；空轴由 sigma 补充 |
| `ssystem` | int | 0 | 0 可分离轴乘积；1 径向。径向模式建议给 slocation；仅轴曲线时实际使用时间表作为径向表 |
| `nlocation` | int[] | 空 | 最多 500 组 `(起始帧, 平面, y, x)`，平铺成数组；每组指定 T 帧连续的 S×S 原始像素区域 |
| `alpha` | float | ftype=0 时 5，否则 7 | 大于 0；采样功率倍数 |

曲线写为 `[频率,值,...]`。每个非空数组须含 0 和 1 端点，频率在 `[0,1]`，值非负，转换为 float32 后位置不得重复；允许乱序，内部排序。被覆盖的曲线也检查原始结构。插值前会按模式变换曲线值，详见[频谱模型](../../knowledge/zh-CN/dfttest/spectral-models.md)。

采样起点满足 `0≤fn≤N-T`，平面必须存在，`0≤x≤plane_width-S`、`0≤y≤plane_height-S`。坐标属于实际平面，不是亮度坐标或块索引；采样没有反射或时间夹取。允许重复组并保留平均权重，也可采样未选中输出的平面。类型 0/1 下非空 nlocation **替换**主 sigma/曲线表；类型 2/3/4 下采样参数仍校验，但不请求采样帧。所有选中平面共享同一采样表。

## 输出与执行

| 参数 | 类型 | 默认 | 取值和作用 |
|---|---|---|---|
| `dither` | int | 0 | 非负；仅选中 UInt8 平面生效。0 普通转换，1 误差扩散，≥2 扩散加确定性噪声 |
| `dither_seed` | int | 0 | 非负；仅 UInt8 且 dither≥2 时参与噪声；0/1 不使用 seed |
| `planes` | int[] | 省略时全部非 Alpha 平面 | 显式 `[]` 不处理任何平面；合法实际索引，重复只选一次 |
| `opt` | int | 0 | 允许 0/1/2/3/8；1 选自有标量核，其余自动 Highway，不强制 ISA，也不强制 FFT 为标量 |
| `threads`、`fft_threads` | int | 省略 | 仅供调用兼容，完全不读取、归一化或保存；不限制为 int32，无执行效果 |
| `fft_backend` | string | 省略 | 仅供调用兼容，内容完全忽略，不选择后端 |

除上述完全忽略的兼容参数外，所有整数及数组整数须在 int32 内，浮点参数须有限且可表示为 float32；不参与当前模式的算法参数仍校验。兼容参数仅受宿主注册类型约束，不进入 DS2 解析或算法配置。固定使用 PocketFFT，滤镜和 FFT 均不自建工作线程；宿主可并发请求帧。没有 DFTTest 跨帧频谱缓存参数。

返回剪辑保持宽高、格式、帧数、帧率；输出帧属性来自源帧 n。未选中平面逐字节复制（AVS 显式模式 1 除外）。`planes=[]` 在配置校验后复制当前帧，不构造模型或 FFT 窗。整数输出最后舍入、裁剪；**float 输出除以 255 后不裁剪**，保留负色度及超出名义范围的值。dither 在所有时空合成完成后只做一次，其他位深忽略合法的 dither 设置。

## 最短示例

```python
import vapoursynth as vs
core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")
clip = core.std.BlankClip(width=128, height=96, length=12,
                          format=vs.YUV420P8, color=[96, 128, 128])
output = core.neo_fft.DFTTest(clip, sigma=8.0, tbsize=3, planes=[0])
output.set_output()
```

以下复用上例源，各自可作为 output：

```python
# 四帧体积，时间重叠两帧。
output = core.neo_fft.DFTTest(clip, tmode=1, tbsize=4, tosize=2)
# 所有频率的共享曲线；仍使用空间默认块和重叠。
output = core.neo_fft.DFTTest(clip, slocation=[0, 4, 1, 12], ssystem=1)
# 从帧 0 开始的三帧、亮度左上 16×16 区域估计噪声。
output = core.neo_fft.DFTTest(clip, nlocation=[0, 0, 0, 0], alpha=5)
```

## 常见错误

减小 sbsize 时同步调整 sosize；例如 smode=1,sbsize=8 不能继续用默认重叠 12。中心模式可以用 sbsize=5，不需改 sosize。T 超过片长、偶数 T 搭配 tmode=0、非法重叠、零能量窗或零 DC 模板会在创建时失败。被消费的 float 样本或中间值出现 NaN/Inf 会导致取帧失败；未选中区域的纯复制不检查样本有限性。

计算流程见 [DFTTest 原理](../../knowledge/zh-CN/dfttest.md)，时间模式见[时间重叠](../../knowledge/zh-CN/dfttest/temporal-ola.md)，量化见 [dither](../../knowledge/zh-CN/dfttest/dither.md)。

## 省略参数与空列表

VS 显式传 `planes=None` 等价于省略参数；AVS 用 `planes=Undefined()`，括号不能省略。两者都不是空列表：`planes=[]` 在两个滤镜中均表示不处理任何平面，校验配置后逐字节复制当前帧。默认处理全部现有非 Alpha 平面；AVS 省略或传 `Undefined()` 时，仍按 `y/u/v/a` 及其默认值选择。显式空列表覆盖 `y/u/v/a`，即使其中含模式 1 或非法模式值。

```python
core.neo_fft.DFTTest(clip, planes=None)  # 使用默认选择
core.neo_fft.DFTTest(clip, planes=[])    # 全部复制
```

```avs
neo_fft_DFTTest(c, planes=Undefined(), y=3, u=2, v=2) # 按 y/u/v/a 选择
neo_fft_DFTTest(c, planes=[])                         # 全部复制
```

## AviSynth

两个滤镜均接受平面 YUVA/RGBA，Alpha 索引为 3。Alpha 为全分辨率，样本缩放按亮度/RGB 处理，不采用色度中心偏移或子采样。未指定选择参数时处理现有非 Alpha 平面并复制 Alpha；`planes=[3]` 可单独处理 Alpha。

AVS 专用整数参数 `y/u/v/a` 按此顺序追加在所有共享参数之后。显式提供 `planes` 时以其为准，忽略这些模式；`planes=[]` 也算显式提供。否则未提供的模式默认 `y=u=v=3, a=2`；1 跳过写入，2 原样复制，3 处理。模式 1 不初始化对应输出平面，其像素内容不作保证，适用于后续只提取有效平面再组合的脚本；它省去该平面的复制，但不缩减输出帧的格式或分配。不存在的平面对应的合法模式无效果；未被 `planes` 覆盖时，其他模式值报错。RGB 中这些名字依次对应 R/G/B/A。`y=2,u=2,v=2,a=2` 在 FFT3D 中也表示全复制。VS 继续使用 `planes` 接口；以独立灰度 clip 保存的 Alpha 可单独送入滤镜。


AviSynth 使用同名参数和顺序；音频与奇偶性来自输入，具体宿主要求见 [API 目录](README.md#avisynth-调用与构建)。

`nlocation`、`slocation`、`ssx`、`ssy`、`sst` 也接受空白、逗号或冒号分隔的数字字符串，可与原生数组混用。`nlocation` 使用十进制 int32 整数；曲线使用十进制浮点数，支持正负号、小数点和科学计数法。小数点固定为 `.`，不受系统区域设置影响。分隔符可混用、连续出现或位于首尾；空串或只含分隔符的字符串等价于空数组。逗号始终分隔数字，不是小数点："1,5" 表示 1 和 5 两个值。方括号、非有限值和带多余字符的数字会被拒绝。解析后仍应用上文的四元组、曲线端点、坐标及取值校验。`planes` 继续使用整数或数组；VS 的这些参数继续使用数字数组。

```avs
LoadPlugin("/path/to/neo-fft.dll")
c=BlankClip(width=128,height=96,length=12,pixel_type="YV12",color_yuv=$608080)
return neo_fft_DFTTest(c,sigma=8.0, tbsize=3, planes=[0])
```

```avs
# 复用上例 c；与 slocation=[0,4,1,12] 相同。
return neo_fft_DFTTest(c,slocation="0:4, 1:1.2e1",ssystem=1)
# 采样也可写为 nlocation="0,0:0,0"。
```
