# BSMPT v3.1.4 与当前主项目差异报告

更新日期：2026-09-04。

## 1. 比较对象与结论摘要

- 旧版：官方 tag `v3.1.4`，commit `0b95ac4dad5bf6eee017f0c7094140b9903054d0`。
- 当前：主项目 HEAD `92f8a9de`，`CMakeLists.txt` 声明 BSMPT 3.3.1；核心官方
  源码基线 `04cb17d1233522f3c423cbd957a8922be037241e`。目录名中的 3.1.8 不是
  实际源码版本。
- 从 v3.1.4 到当前官方基线共 37 个提交；核心范围 69 个文件，约
  `+8000/-2298` 行。当前 HEAD 相对 `04cb17d1` 没有额外主项目核心 C++ 差异。
- 理论框架没有更换，但有限温导数、AE Daisy 导数、声速、相追踪和温度求解包含
  实质修正，不能把两个版本预期为逐位一致。
- 八个已完成 R2HDM 分层点中，状态与 transition history 全部一致；强信号点总
  SNR 差 5.82%，两个极弱信号点相对差 89.58% 和 161.78%，但绝对 SNR 仅
  `10^-26` 和 `10^-15`。

## 2. 构建与 double-free

v3.1.4 的 `conanfile.py` 明确指出：Linux 下 BSMPT 与 libcmaes 的 vectorization
配置不一致会触发 double free。第一次仅传
`-DBSMPTUseVectorization=OFF`，当前 Conan toolchain 仍注入 `-ftree-vectorize`，
高 SNR 点在初始化后复现 `double free or corruption`，且没有生成 TSV。

正式旧版对照树 `build-v3.1.4-novec/` 额外使用：

```text
-fno-tree-vectorize -DEIGEN_DONT_VECTORIZE
```

修正后同一点完整运行 240.667 s，不再崩溃。两次构建均只存在于本比较目录，未修改
v3.1.4 源码、当前主项目或 `install/bin/CalcGW`。

## 3. 理论与物理实现变化

### 3.1 未改变的框架

两个版本仍使用 Coleman-Weinberg 零温一圈势、玻色/费米有限温热函数，以及
Parwani 或 Arnold-Espinosa (AE) 热质量方案。碰撞、声波、湍流三类 GW 谱的总体
框架和输入量仍是 `T*`、`alpha`、`beta/H`、声速、壁速和 `g*`。

### 3.2 热函数数值实现

提交 `0564f24e` 将 `Class_Potential_Origin::boson/fermion` 导数路径的
`JbosonNumericalIntegration/JfermionNumericalIntegration` 换成分段插值/渐近展开。
这是数值实现改变，不是热势定义更换；主要影响速度、连续性和导数噪声，偏差方向
不能普遍预言。

### 3.3 显式温度导数修正

v3.1.4 的热函数温度导数相当于假设质量平方不显含温度。`0564f24e` 增加
`dm²/dT` 项，把

```text
4 T^3 J - 2 T m² J'
```

修正为

```text
4 T^3 J - (2 T m² - T² dm²/dT) J' .
```

它会直接改变 `dVeff/dT`、能量差和 `alpha`；符号取决于具体质量谱。

### 3.4 AE Daisy 与规范玻色子导数

`0564f24e`、`431df3b6` 补全 Daisy 的场导数和温度导数，并修正 Parwani 规范玻色子
导数误用 Higgs 常数的问题。这些修正会改变极小值、势垒、alpha 和 beta/H，尤其
可能放大边缘弱信号的相对差异。

### 3.5 声速与壁速

`431df3b6` 修改 `BounceSolution::CalculateWallVelocity()`：旧版依赖相邻温度有限
差分，当前使用 `VEff(...,-1)` 的温度导数并通过 `NablaNumerical` 求二阶导数。
因此 `cs_f/cs_t`、壁速、声波效率和最终 SNR 都可能变化，方向依参数点而定。

### 3.6 相变温度与终止性

- `0a28d291` 修复 `TrackPhase` 的 `dT==0` 无限循环、终点步长、`Tc==T_low` 和多步
  递归条件，并将 multistep mode 改为强类型枚举。
- `48781cff` 重做高温全局真空选择和共存相区域组织，可能改变相编号、共存对与 Tc。
- `906b818e` 将 false-vacuum fraction 温度搜索改为最多 100 次二分，并在 action
  重采样后重算渗流/完成温度。边缘点可能从继续迭代变为 `NotMet/-1`。

上述属于稳定性/算法修正，而非新相变理论。

### 3.7 其他重要变化

- `65b5e1c2` 按积分方向设置 GSL 初始步长，避免 ODE 方向错误。
- `9aebad88` 将近 Tc action 缓存阈值从 `1e-3` 收紧至 `1e-4` 并改善日志。
- `01578231` 统一模型旋转矩阵符号约定；`0cba998a` 调整小元素容差。
- `2a41b3d1` 修复 Z2 对称性识别，可能影响对称性约化和 VEV 扇区。
- `41a98730` 新增 RxSM；对本报告的 R2HDM 点无直接作用。

## 4. 八点数值结果

逐点原始统计见 `version_result_summary.tsv`。相对偏差统一定义为
`(current-v3.1.4)/abs(v3.1.4)`。

| 样本 | GW 状态 | history | alpha | beta/H | 总 SNR | runtime 3.1.4→当前 |
|---|---|---|---:|---:|---:|---:|
| broad A | success→success | 相同 | +25.14% | -0.0023% | +161.78% | 173.924→175.116 s |
| broad B | nan→nan | 相同 | 不适用 | 不适用 | 不适用 | 17.371→17.311 s |
| broad C | nan→nan | 相同 | 不适用 | 不适用 | 不适用 | 9.190→9.344 s |
| NLO valid | success→success | 相同 | +15.64% | -0.134% | +89.58% | 78.702→136.795 s |
| NLO invalid | nan→nan | 相同 | 不适用 | 不适用 | 不适用 | 0.336→0.312 s |
| high SNR type-3 | success→success | 相同 | +3.12% | -1.184% | +5.824% | 240.667→238.323 s |
| Yukawa type-2 | success→success | 相同 | -25.85% | -0.0256% | -57.46% | 220.381→220.126 s |
| Yukawa type-4 | success→success | 相同 | -25.85% | +0.0192% | -57.93% | 223.005→222.972 s |

状态/history 翻转为 `0/8`。五个 success 点均保持成功，三个失败/拒绝点保持相同
失败状态。高 SNR 点还表现为：`Tstar +0.0556%`、`Treh +0.4467%`、`cs_f +0.239%`；
温度变化小于 0.5%，主要偏差来自热导数、alpha 和声速链。

broad A 与 NLO-valid 的总 SNR 分别从 `8.215e-16→2.151e-15` 和
`2.561e-26→4.855e-26`。相对差很大，但绝对值极弱，不应与可观测强信号点的
5.82% 混为一谈。

八点 runtime 合计为 963.576 s→1020.299 s，当前版本在该小样本慢 5.89%。这不是
通用性能结论：样本小、路径构成差异大，而且 v3.1.4 必须显式关闭 vectorization
才能与当前 libcmaes 安全配合。高 SNR 单点当前反而快 0.97%。

## 5. 结论与适用边界

当前版本不是单纯优化版，而是包含物理导数和相追踪修正的后续官方版本。对稳定
高 SNR 示例，状态不变且 SNR 偏差约 5.82%；对弱信号边缘，alpha/SNR 相对偏差可
达几十至上百个百分点，即使状态仍相同。因此不能用 v3.1.4 的 SNR 对当前版本做
统一常数缩放，也不能假定偏差总是同方向。

本轮八点没有状态翻转，但这不是全参数空间保证。最容易出现版本级状态差异的区域
仍是高温真空选择、共存相边界、`Tc==Tlow`、false-vacuum fraction 未收敛和 Z2
约化边界。扫描迁移时应以当前版本重新计算这些点，而不是复用 v3.1.4 的状态。
