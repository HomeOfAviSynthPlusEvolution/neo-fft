# UInt8 量化与确定性 dither

[目录](../README.md)

DFTTest 仅对选中 UInt8 平面启用 dither。它位于完整时空合成、裁取之后，不改变频谱功率或窗。高位深整数和 float 对所有合法模式仍走普通转换；负模式或负 seed 在不生效时也报错。

[API](../../../api/zh-CN/dfttest.md) · [样本尺度](../shared/sample-domain.md)

## 扫描与误差

每帧、每平面先清零两行误差。每一行都从左到右，行从上到下；不是蛇形，不在块边界重置。设 E 为重建后的 8-bit 幅度，e 为当前传入误差：

```text
dither=0: v = E+0.5
dither=1: v = (E+e)+0.5
dither>=2:
    scale = float(dither-1)+0.5
    off = scale*0.5
    v = (((E+u*scale)-off)+e)+0.5
D = clip(trunc(v),0,255)
error = E-D
```

error 用原 E 减去已裁剪输出 D，**不是**加过噪声和传入误差的 v−D。依次把误差送往下一行左、下一行同列、当前行右、下一行右，权重为 3/16、5/16、7/16、1/16。图外贡献丢弃、不重新归一化；每帧独立。

例一行两个像素 E=[10.4,10.4]、dither=1：首像素 D=10，误差 0.4，给右邻 0.175；第二像素 v=10.4+0.175+0.5=11.075，D=11。第二像素自身误差为 −0.6，不是 0.075。

## 噪声由坐标确定

n 是输出帧号，p 是实际平面索引，x/y 是该可见平面坐标；与 planes 数组顺序、块号或线程身份无关。

```text
mix(v):                              # uint32, modulo 2^32
    v = (v ^ (v >> 16)) * 0x7feb352d
    v = (v ^ (v >> 15)) * 0x846ca68b
    return v ^ (v >> 16)
h = mix(uint32(seed) ^ 0xa511e9b3)
for coordinate in [n,p,y,x]: h = mix(h ^ uint32(coordinate))
u = float(h >> 8) * 2^-24
```


seed=n=p=y=x=0 时 h=0x642a1d14，h>>8=6564381，u=6564381/16777216。dither=2 时 scale=1.5,off=0.75，所以该像素的扰动为 `1.5u−0.75`，约 −0.1631。相同坐标的 u 不因先请求了另一帧而改变。

模式 0/1 不使用 seed，也不推进隐藏随机状态。≥2 的噪声定义是本项目确定性契约，不承诺匹配旧实现的工作区 RNG。由于浮点运算顺序会影响半 LSB 分支，不能把 `E+u*scale-off` 重排为 `E+(u*scale-off)`。噪声可 SIMD 生成，但扩散扫描有依赖，当前在调用线程内顺序执行。

有限性检查与安全饱和先于危险整数转换；失败不发布半张帧。固定算术配置下，重复、逆序和并发请求应得到相同输出，参见[执行精度](../shared/execution-precision.md)。
