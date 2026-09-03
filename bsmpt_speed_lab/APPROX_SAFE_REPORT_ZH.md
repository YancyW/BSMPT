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
