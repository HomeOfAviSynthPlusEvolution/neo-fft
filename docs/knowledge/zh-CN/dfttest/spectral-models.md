# DFTTest 频谱增益、曲线与采样

[目录](../README.md)

令 X 为已加窗的时空频谱，G=FFT(255h)。开启 zmean 时，`g=Re(X[0])/Re(G[0])`、`M=gG`、`R=X−M`；否则 M=0,R=X。所有滤波作用于残差 R，最后加回 M。

[DFTTest 主文](../dfttest.md) · [完整参数](../../../api/zh-CN/dfttest.md)

## 五种频谱规则

设 `p=Re(R)²+Im(R)²`、ε=1e−15，wscale=1/Σh²。主模型值为 v：类型 0/1 的 A=v/wscale，类型 2/3/4 的 A=v。B 对 sigma2 用同样除数；L=pmin/wscale、H=pmax/wscale。

| ftype | 残差增益 |
|---|---|
| 0 | a=max((p−A)/(p+ε),0)，gain=a^f0beta |
| 1 | p<A 时 0，否则 1；等号保留 |
| 2 | A |
| 3 | L≤p≤H 时 A，否则 B；两端均包含 |
| 4 | q=p+ε，gain=A√(qH/((q+L)(q+H))) |

类型 0 若 float 距离 `abs(f0beta−1)<0.00005`，直接用 a；否则若 `abs(f0beta−0.5)<0.00005` 用 √a；其余用 pow。比较严格小于。公式中 epsilon 的位置与 FFT3D 不同，不应共用一个经过“化简”的 Wiener 式。

例 R=3+4i、p=25、A=9：类型 0、指数 1 得 gain≈0.64，结果≈1.92+2.56i；指数 0.5 得 gain≈0.8。类型 1 在 A=25 时保留 R，A>25 才清零。类型 3 的 L=25,H=100,A=2,B=0.5 在 p=25 时用 2。

## 曲线如何成为 v

所有曲线为空时直接使用 sigma。否则 d=(T>1?1:0)+(B>1?2:0)。长度 L>1 的轴频率为 `f(i,L)=min(i,L−i)/floor(L/2)`，长度 1 取 0；T=3 映射为 `[0,1,1]`。

先将曲线值取 e 次幂，**然后**线性插值：共享 slocation 在 ssystem=0 取 e=1/d，在 ssystem=1 取 e=1；仅轴曲线时两个系统都取 e=1/d。缺省轴曲线用 sigma^e。精确 knot 直接取该值。

ssystem=0 用时间、Y、X 三个插值表相乘，单元素轴贡献 1。ssystem=1 算 `rho=sqrt((ft²+fy²+fx²)/d)`，使用时间表在 rho 的值；仅轴曲线时 ssx/ssy 不提供径向表，所以一般径向需求应给 slocation。d=0 时取 slocation 的 DC，否则 sst 的 DC，否则 sigma，指数为 1。

例如 T=1,B=3，d=2，共享可分离曲线 `[0,4,1,16]` 先变为 `[0,2,1,4]`。在 fx=1、fy=0 的 bin，v=4×2=8；不是先插值原始曲线后随意开方。径向共享模式则在 rho=1/√2 直接插值 4 到 16。

## 采样表替换主模型

对类型 0/1，非空 nlocation 优先于全部主 sigma/曲线。采样窗 h2 使用原始空间、时间窗，两轴都不做重叠归一化，但保留 1/√(TB²)。每个 tuple 消费连续 T 帧的原始矩形，按正常 DFTTest 幅度单位变换；开启 zmean 时用 G2=FFT(255h2) 移除采样窗形均值。

```text
wscale2 = 1/sum_float(h2*h2)
Psum[k] = sum_in_tuple_order(abs(R_sample[k])**2)
A[k] = Psum[k] * ((1/M)*(wscale2/wscale)*alpha)
```

M 包含重复 tuple。A 已经是校准后的功率，不能再次除 wscale、平方或乘 T。所有输出平面共享它；采样不改变 sigma2/pmin/pmax。例两份残差功率 4、8，M=2，wscale2/wscale=0.5，alpha=5，则 A=15。

原始数组即使被覆盖也校验；被覆盖的派生曲线不必计算。负值、重复 float32 knot、缺失端点、越界采样或非有限活动功率会报错。采样均值与普通滤波均值各用自己的模板，不能混用。所得结果遵守[执行精度](../shared/execution-precision.md)。
