# 样本尺度与输出转换

[目录](../README.md)

滤镜都用 float32 工作缓冲，但像素数值的含义不同。这里的 bits 指整数输入位深，q 指进入分析窗前的幅度；没有电视范围到全范围的自动映射。

| 输入 | FFT3D 的 q | DFTTest 的 q |
|---|---|---|
| 整数亮度、RGB、Alpha | sample | sample / 2^(bits−8) |
| 整数 YUV 色度 | sample − 2^(bits−1) | sample / 2^(bits−8)，不减中点 |
| float32，所有平面 | sample | sample × 255 |

例如同一整数亮度 128、512、32768（8/10/16 位）在 DFTTest 都变成 128；FFT3D 保持各自尺度，并将 sigma 相应乘 1、4、256。float 亮度 0.5 在 DFTTest 变成 127.5，在 FFT3D 保持 0.5，其 sigma 除以 255。10-bit 中性色度 512 在 FFT3D 为 0，在 DFTTest 为 128。float 色度 −0.25 在 DFTTest 为 −63.75；两者都不再减色度中点。

FFT3D 的 sigma 是标准差，形成噪声功率时平方；DFTTest 的 sigma 是类型相关功率或增益，不能平方，也不能跨滤镜用同一数值解释。

## 输出在哪里裁剪

最终合成值记为 z。FFT3D 整数依次计算 `(z+0.5)+base`，安全截断并裁到 `[0,2^bits−1]`；base 仅在整数色度时为中点。float 裁到 `[0,1]`，包括 float 色度。因此负色度经过 FFT3D 处理可变成 0。

DFTTest UInt8 普通转换为 `clip(trunc(z+0.5),0,255)`；高位深为 `clip(trunc(z*2^(bits−8)+0.5),0,2^bits−1)`；float 为 `z/255`，不裁剪。UInt8 的可选 [dither](../dfttest/dither.md)替换最后的普通转换。

例如 DFTTest 合成 z=128.25，普通 u8 输出 128，10-bit 输出 513，float 约 0.502941。这说明不要先量化到 u8 再恢复高位深。饱和发生在任何危险的整数转换之前；没有对块中间结果逐块裁剪。

## 复制、格式和有限性

两个入口支持固定平面 GRAY/YUV/RGB 的 1 或 3 平面，AviSynth 另支持四平面 YUVA/RGBA；整数 8/10/12/14/16 或 float32。重复 planes 索引只选一次。两者默认选择所有现有非 Alpha 平面；显式 `planes=[]` 均选无，校验配置后全部复制。VS `planes=None` 与 AVS `planes=Undefined()` 表示省略，不等于空列表；AVS 省略时仍使用 y/u/v/a 及其默认值。AVS 显式索引 3 选择 Alpha，其尺寸和样本缩放按亮度/RGB 处理，不采用色度子采样或中心偏移。AVS y/u/v/a 仅在未显式提供 planes 时生效，默认 3/3/3/2，1 跳过写入、2 复制、3 处理。模式 1 不构造该平面的处理计划，并在普通输出和 Kalman 输出（含帧 0）中都跳过写入，不保证输出像素内容；完整输出帧仍会分配。其余未选平面、选中平面的 FFT3D ROI 外部以及 Kalman 输出帧 0 逐字节复制，保留原样本位模式。

活动计算所消费的浮点样本与中间值要求有限。DFTTest 的采样模型可以消费未选中输出平面的指定 patch，因此该 patch 也须有限。复制区域不因包含 NaN 而自动失败。整数与 float 输出不改变尺寸、帧数或帧率，属性取自源帧 n。

参见 [FFT3D API](../../../api/zh-CN/fft3d.md)、[DFTTest API](../../../api/zh-CN/dfttest.md)和[执行精度](execution-precision.md)。
