# FFT3D 噪声来源与增强次序

[目录](../README.md)

[FFT3D 主文](../fft3d.md)给出空间／时间 Wiener；这里说明它消费的噪声功率从哪里来，以及增强为什么不能任意交换顺序。参数范围见 [API](../../../api/zh-CN/fft3d.md)。W=bw,H=bh；所有频点公式只遍历实际半谱。

## 三种噪声来源

创建时选择：pfactor>0 用采样；否则四个有效 sigma 不全相等时用频率曲线；否则均匀模型。输入 sigmas 先按整数 `2^(bits−8)` 或 float `1/255` 缩放。均匀空间功率为 `P=W*H*sigma_eff²`，不使用窗能量代替块面积。

频率曲线如下，式中的 sigma 均已缩放：

```text
norm = 1/(W*H)
fy = (H - 2*abs(y-floor(H/2)))/H
fx = x/(floor(W/2)+1)
f = sqrt((fx*fx+fy*fy)/2)
a = sqrt(0.5)/4; b = sqrt(0.5)/2
s = sigma4 + (sigma3-sigma4)*f/a             if f<a
    sigma3 + (sigma2-sigma3)*(f-a)/(b-a)    if a<=f<b
    sigma + (sigma2-sigma)*(1-f)/(1-b)      otherwise
P = s*s/norm
```


不要将 fx 分母改成 floor(W/2)，或将 fy 换成常见的对称频率定义。奇数 H 的行为尤其不同。普通 T 帧 Wiener 使用 T*P；保存的基础 P 不随第一次请求的首尾回退而改变。

## 采样、评分和预览

从夹取后的 pframe 取正常反射／分析窗块谱 X，用 `g=degrid*Re(X[0])/Re(G[0])` 去除窗形网格，degrid=0 时直接取 R=X。设 `fy=2min(y,H−y)/H`、`fx=2x/W`、`w=(fy²+fx²)/(fy²+fx²+pcutoff²)`，评分为 Σ|R|²w，不乘半谱重数。

手动 px/py 是块索引；两者同时为 0 则在 bx=2..Nx−3、by=2..Ny−3 搜索，按 Y/X 顺序保留第一个最小评分。选定后普通噪声表为 `P=pfactor*|R|²*w`。例如 |R|²=100,w=0.4,pfactor=2 得 P=80；T=3 使用 240。Kalman 是例外：正 pfactor 仅选择采样，R 使用未乘 pfactor 的基础功率，即本例 40。

采样或非均匀模型开启时，pshow 可以覆盖其他处理路径。它从当前 n 用正常窗选择块，再以矩形分析／合成窗仅重建该块，不降噪、不增强、不读 pframe，不发布采样模型。ROI 中其余块贡献为零；整数色度恢复中点。不要把预览当作降噪后的噪声图。

## 增强的频率权重与功率增益

对待增强空间谱 Z 重新求 g，R=Z−gG。smin/smax 与 sigma 同样按格式缩放，ht **不缩放**。使用：

```text
dy = y if y<floor(H/2) else H-y
d2 = dy*dy*svr*svr/floor(H/2)^2 + x*x/floor(W/2)^2
Ws = 1-exp(-d2/(2*scutoff*scutoff))
rawWh = exp(-0.7*d2*hr*hr)-exp(-d2*hr*hr)
Wh = rawWh/max(rawWh)
A = smin_eff^2/norm; B = smax_eff^2/norm; C = ht^2/norm
q = abs(R)^2+1e-15
S(q) = 1 + sharpen*Ws*sqrt(q*B/((q+A)*(q+B)))
D(q) = (q+C)/((q+C)+dehalo*Wh*q)
output = R*S(q)*D(q)+gG
```


强度为零时对应增益精确旁路，dehalo=0 不计算未使用的 Wh 归一化；活跃 Wh 必须有正且有限最大值。B=0 时锐化增益为 1。H=5 的 dy 为 `[0,1,3,2,1]`，不是通常的对称最小距离。

两增益组合使用同一残差功率。例如 q=100,A=16,B=400,C=2500,Ws=Wh=1,sharpen=0.5,dehalo=0.2 时 S≈1.41523,D≈0.992366，总幅度因子约 1.40442。这个单 bin 例子把权重视为已给定；不表示每个实际频点都能同时取得 Ws=Wh=1。

## 哪一步使用这些增益

| 路径 | 次序 |
|---|---|
| 有效单帧、均匀噪声 | 同一原始 R/q 下融合 Wiener×S×D，再恢复网格 |
| 有效单帧、曲线或采样 | 先 Wiener 和网格恢复，再对结果重新求增强 g/q |
| T=2..5 | 时间滤波后恢复当前空间谱，再增强 |
| bt=−1 | 原始空间谱直接增强，无 Wiener |
| bt=0 | 在 L 的输出副本上增强，不反馈递推 |
| 有效预览 | 不使用这些增益 |

假设 Wiener 使残差幅度减半，功率会变成原来的四分之一。先滤波再增强所使用的 S(q/4) 一般不同于融合路径的 S(q)，所以两条路径不能仅因噪声表碰巧恒定就合并。各分支固定其运算次序；参见[执行精度](../shared/execution-precision.md)和 [Kalman](kalman.md)。
