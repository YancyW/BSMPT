# CalcGW 非严格安全回退版研究记录

## 目标与边界

该模式用于加速参数扫描，不替代严格版的最终发表结果。候选只增加
`BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0=1`：R2HDM 的 `diff==0` lepton 9x9
矩阵按三个解析块求谱，其余均继承已经过零容差验证的 exact-fast 配置。解析实现
会检查完整稀疏结构和对称性，结构不符时自动回退完整 Eigen。

验收目标是：已接受的近似结果中，所有状态与 transition history 不变，最终及各
分量 SNR 相对偏差不超过 10%。由于相变搜索会放大机器舍入，有限样本测试不能对
整个连续参数空间给出数学保证，因此正式入口采用“近似首遍 + 风险点严格回退”。

## 直接近似测试

| 数据组 | 行数 | 状态/历史 | 最大 SNR 相对偏差 | exact-fast / s | approximate / s |
|---|---:|---|---:|---:|---:|
| NLO 边界内侧 | 5 | 全部相同 | 1.016% | 183.053 | 154.598 |
| NLO 边界外侧 | 4 | 全部相同 | 0 | 约 0.7 | 0.725 |
| 广域 A | 10 | 1 行发生 bounce 状态翻转 | 不作为接受组 | 387.929 | 327.026 |
| 广域 B | 10 | 全部相同 | 11.60%（总 SNR 11.58%） | 825.687 | 737.946 |
| 广域 C | 10 | 全部相同 | 1.107% | 556.091 | 508.728 |
| 高 SNR × Yukawa 2/3/4 | 3 | 全部相同 | 4.360% | 220.162 | 171.789 |

高 SNR type 3 的总 SNR 为 151.086→151.741，相对偏差 0.433%。该三行批次
耗时降低 21.97%。C 组包含约 313 s 长尾，整体只降低 8.52%，说明近似收益取决于
lepton 谱在总耗时中的比例。

两个直接反例证明不能盲目启用解析谱：A 组第九数据行从 bounce `success`、后续
`not_met` 变为 bounce `failure`、后续 `nan`；B 组第四数据行的总 SNR 从
`2.929e-27` 变为 `3.269e-27`，相对偏差 11.58%。

## guarded 策略

单点入口 `run_calcgw_approx_safe.sh` 先调用 `run_calcgw_approx_guarded.sh`，然后由
`approx_needs_fallback.py` 检查输出。以下任一条件触发 exact-fast 重算并覆盖输出：

- 任意状态为 `failure`、`nan`、`not_met` 或 `not_set`；
- 总 SNR 缺失、非有限或绝对值低于 `BSMPT_APPROX_SNR_FLOOR`（默认 `1e-20`）；
- 总 SNR 落入用户关心 cut 的安全带。cut 用逗号分隔的
  `BSMPT_APPROX_SNR_CUTS` 指定，半宽由 `BSMPT_APPROX_CUT_MARGIN` 控制，默认
  15%。

默认地板在当前 42 行对照中接受 6 个稳定 GW 点、回退 36 个风险或极弱信号点；
被接受点的最大 SNR 分量相对偏差为 4.36%，没有超过 10%。已知 A 状态翻转和 B
弱信号超限点均被检测器捕获。端到端测试中，A 状态翻转点自动严格重算后与严格
基准所有非 runtime 字段完全一致；高 SNR type 3 则接受近似结果，耗时约 62.0 s，
总 SNR 偏差 0.433%。

在近似首遍之前，入口还用 `approx_input_prefilter.py` 查询只会增加严格计算的已知
风险注册表。当前仅登记 A 组精确反例和经过四点 `±1e-5` 扰动验证的 B 组危险
邻域；命中时直接运行 exact-fast，未命中或参数无法解析时保持原有 guarded 流程。
42 点离线回归只命中预期 A/B 两点，B 的四点扰动云 4/4 命中。A 点真实集成测试
确认日志为 `direct exact`，最终非 runtime 字段与严格参考逐位一致，并避免了约
41.337 s 的必然回退首遍；B 锚点可避免约 101.822 s。注册表不是参数空间分类器，
不得把未经 exact/approx 邻域验证的新范围加入其中。

## 使用方式

该入口要求一次只计算一个输入行，以便风险点独立回退：

```bash
bsmpt_speed_lab/run_calcgw_approx_safe.sh R2HDM input.tsv output.tsv 2 2
```

如果扫描以 SNR=10 和 100 为分类边界，可使用：

```bash
BSMPT_APPROX_SNR_CUTS=10,100 \
  bsmpt_speed_lab/run_calcgw_approx_safe.sh R2HDM input.tsv output.tsv 2 2
```

风险回退会先付出一次近似计算再运行严格计算，因此它在大量失败/弱信号点上可能
比直接 exact-fast 更慢。它适合预期存在较多明确强信号、且需要保护 cut 附近点的
扫描。最终入选点仍建议用 `run_calcgw_exact_fast.sh` 复核。

## central2 guarded 首遍

进一步测试把 `BSMPT_USE_CENTRAL2_GRADIENT=1` 与解析 lepton 谱组合，仅用于
guarded 首遍。42 行完整矩阵（NLO 内/外侧、A/B/C、高 SNR Yukawa）中，原始
首遍合计从 exact-fast 的 2173.759 s 降至 1597.228 s，减少 26.52%。默认 guard
接受 6 行、严格回退 36 行；六个接受点的最大 SNR 分量相对偏差为 3.89%。

central2 产生了一个额外的重要反例：B 组第八数据行从严格版 bounce `failure`、
GW=`not_set` 变为近似版完整 GW=`success`。该近似信号只有 `3.21e-28`，因此被
默认 `1e-20` SNR 地板捕获并强制严格重算。A 组原有的反向状态翻转也仍被风险
状态检测捕获。C 组 10/10 状态与历史不变，最大 SNR 偏差 0.328%，耗时降低
34.84%；高 SNR 三型最大偏差 3.89%，耗时降低 29.34%。

基于这些结果，`run_calcgw_approx_safe.sh` 默认选择
`run_calcgw_approx_central2.sh` 作为首遍。可通过 `BSMPT_APPROX_FIRST_PASS` 指回
同目录内的 `run_calcgw_approx_guarded.sh`，恢复 analytic-only 首遍。直接调用
central2 wrapper 仍不安全，也不属于可接受结果入口。

### 假阳性边界邻域

围绕 B 组 central2 假 GW success 点，对除 Yukawa type 外的七个输入同时施加
`-1e-5、-3e-6、+3e-6、+1e-5` 相对扰动。四点 exact-fast/central2 结果为：

- 三点 exact-fast 为 bounce `failure`、GW=`not_set`，central2 均误报完整
  GW=`success`；近似总 SNR 为 `3.21e-28--1.15e-27`，全部低于默认地板并回退；
- 一点 exact-fast 为完整 GW=`success`，central2 的 GW 为 `failure/nan`，由风险
  状态规则回退；
- 因此该邻域 4/4 错误近似结果均不会成为 guarded 最终输出。

这一邻域的原始 central2 总耗时为 663.720 s，exact-fast 为 689.692 s，只减少
3.77%，前三点 central2 甚至更慢。触发回退后总成本必然高于直接 exact-fast。
因此安全近似版并不承诺在边缘区加速；它以边缘区严格回退换取稳定强信号区约
29% 的已测首遍收益。邻域输入和输出以
`central2_false_positive_neighborhood_4*` 命名。

## central2 + adaptive threshold

在 central2 guarded 首遍上增加 64 点 adaptive exact-solution threshold 搜索。
同一套 42 行矩阵中，原始首遍合计为 1106.233 s，对应 exact-fast 2173.759 s，
减少 49.11%。默认 guard 仍接受相同 6 行、回退 36 行；接受点的最大 SNR 分量
偏差为 3.72%，未超过 10%。

- NLO 边界内侧 5/5 状态与历史不变，最大 SNR 偏差 1.95%，耗时
  183.053→82.520 s；外侧 4/4 保持拒绝；
- 广域 A 保留同一个已知风险翻转并被回退，耗时 387.929→196.043 s；三个接受点
  最大 SNR 偏差 0.0591%；
- 广域 B 保留同一个低 SNR 假阳性并被回退，耗时 825.687→466.412 s；
- 广域 C 10/10 状态与历史不变，最大 SNR 偏差 0.328%，耗时
  556.091→273.268 s；
- 高 SNR Yukawa 2/3/4 三行状态与历史不变，最大 SNR 分量偏差 3.72%，耗时
  220.162→87.213 s；
- `multistepmode=0/1/2/auto` 均保持状态与历史，最大 SNR 偏差 3.68%；模式 0、1、
  2、auto 分别为 27.163、27.254、89.256、27.487 s，对应 exact-fast
  70.471、72.175、143.627、71.760 s。

在四点假阳性扰动邻域中，adaptive 首遍仍产生三次低 SNR 假 success 和一次
failure，但 4/4 全被 guard 回退；首遍合计 479.302 vs 689.692 s。基于上述门槛，
`run_calcgw_approx_safe.sh` 的默认首遍升级为
`run_calcgw_approx_central2_adaptive.sh`。analytic-only 和 central2-only wrapper
继续保留用于消融；严格 wrapper、严格二进制和项目主程序不变。

## adaptive + raster500

将 bounce `dV/dl` raster 从 1000 降为 500 后，同一 42 行矩阵首遍为
940.724 s，对应 exact-fast 2173.759 s，减少 56.72%；相比未缩减 raster 的
adaptive 首遍 1106.233 s 再快 14.96%。guard 仍接受相同 6 行、回退 36 行，
接受点最大 SNR 分量偏差为 3.89%。

NLO 内侧五点最大 SNR 偏差 5.05%，耗时 82.520→67.933 s；外侧四点保持拒绝。
A/B/C 分别为 163.125、400.762、234.463 s，对应 adaptive 的 196.043、
466.412、273.268 s；A/B 仍只有已知且可捕获的状态翻转，C 组状态与历史不变。
高 SNR 三型为 73.695 s，对应 exact-fast 220.162 s，最大 SNR 偏差 3.89%。

`multistepmode=0/1/2/auto` 分别为 22.505、22.682、80.338、22.310 s，四种模式
状态与历史均不变，最大 SNR 偏差 3.68%。因此 safe runner 默认首遍进一步升级为
`run_calcgw_approx_c2_adaptive_raster500.sh`。更保守的 wrapper 均保留用于回退与
消融，严格路径不变。

### raster250 / raster100 下探

raster250 在高 SNR 和危险边界上相对 raster500 只再快约 6%，NLO 内侧五点
67.933→64.824 s（4.58%）；高 SNR 最大 SNR 偏差 3.87%，NLO 最大偏差 4.78%。
raster100 相对 raster250 仅再快约 3%：高 SNR 22.842→22.028 s，NLO 五点
64.824→62.845 s；最大 SNR 偏差分别 3.81% 和 5.10%。两者的边际收益不足以
补偿更低插值密度带来的未覆盖风险，因此不升级默认，保留为实验 wrapper。

### 当前近似热点

在 raster500 高 SNR type-3 点开启计数和 VEff timing，约 200 万次 VEff 调用中，
累计主要耗时为 quark masses 8.68 s、Higgs masses 6.57 s、fermion thermal
1.98 s、CounterTerm 1.85 s、boson thermal 1.59 s、gauge masses 1.49 s；解析
lepton 已降至 0.64 s。历史 quark 6x6、RK5 fixed/workspace、directional dV/dl
候选均更慢或改变边界，解析 gauge 又被当前 paired gauge 路径旁路，故不重新加入。

### analytic gradient 与 adaptive32 消融

在 adaptive64+raster500 基础上，分别测试了解析梯度替代 central2，以及将 adaptive
grid 从 64 减为 32。解析梯度的强信号 type-3 点状态与跃迁史不变，总 SNR 偏差
0.642%，耗时 26.709 s；但它慢于当前 central2 的 24.326 s。A 组风险点仍发生
bounce 状态差异并被 guard 捕获；B 组危险点与严格版同为 failure。更重要的是，
NLO 内侧五点虽然状态与历史一致，但最大 SNR 分量偏差达到 13.19%，超过 10% 门槛，
且五点均因极低 SNR 回退；其 86.047 s 首遍也慢于当前 67.933 s。因此解析梯度
候选不升级默认。

adaptive32 在强信号点和 NLO 五点的所有非 runtime 输出与 adaptive64 完全相同，
但实测分别为 26.314 vs 24.326 s、96.662 vs 67.933 s，没有速度收益，也不升级
默认。两项消融的 wrapper 与原始 TSV 留在实验目录，便于复查；严格路径未变。

### raster400 中间点

raster400 用于检查 raster500 与 raster250 之间是否存在更好的折中。NLO 内侧五点
状态/历史一致、最大 SNR 分量偏差 4.77%，耗时 66.377 vs raster500 的 67.933 s；
高 SNR 三型状态/历史一致、最大偏差 3.84%，但总耗时 73.367 vs 73.695 s，仅快
0.45%。A、B 两个已知危险点仍分别由失败状态和低 SNR 地板捕获。完整 C 组 10/10
状态/历史一致、最大 SNR 偏差 0.352%，但耗时 238.339 s，反而慢于 raster500 的
234.463 s。收益不可重复，因此 raster400 也不升级默认。

### PGO 复用审计

现存 `build-pgo-gcc` 仍是 `-fprofile-generate` instrumentation 构建，并非可测速的
profile-use 成品。虽然保留 35 个 `.gcda`，但当前源码时间戳晚于训练 binary，无法
证明画像匹配。因此没有运行该 binary，也不把历史 2–5% PGO 收益叠加到安全版；
未来若重启，必须从当前源码和广域训练集完整重做。
