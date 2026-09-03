# CalcGW 优化研究逐轮日志

本文件是 `bsmpt_speed_lab/` 的反查入口。每轮研究无论成功、失败或仅得到负结果，
都必须记录：日期、思路、相对基线、样本、判定标准、结果、结论、产物以及 Git
提交/推送状态。详细数值解释见 `APPROX_SAFE_REPORT_ZH.md`，严格优化历史见
`OPTIMIZATION_REPORT_ZH.md`。

## 固定边界与验收规则

- 所有修改、wrapper、构建和输出只能位于 `bsmpt_speed_lab/`。
- 项目主程序、严格分支/标签和严格二进制不得修改或重建。
- 同时最多运行两个 CalcGW，避免超过 WSL 24 GB 内存限制。
- guarded 接受结果必须保持全部 status 与 `transition_history`，所有 SNR 分量相对
  偏差不超过 10%；失败、非有限值、极弱信号和 cut 安全边缘自动 exact-fast 回退。
- 被淘汰候选也保存 wrapper、原始 TSV 与结论，避免重复试验。
- 每轮完成后更新本日志；形成可复用结论时提交并推送特殊研究分支。

## 提交状态总表

| 轮次 | 内容 | 结论 | commit | GitHub 状态 |
|---|---|---|---|---|
| S0 | exact-fast 严格快照 | 接受并冻结 | `f892f02e` | 已推送严格分支并创建标签 |
| A1 | analytic-lepton guarded 基线 | 保留为消融 | `f4ffda0e` | 已推送 |
| A2 | central2 guarded | 当时升级默认 | `95ff109e` | 已推送 |
| A3 | adaptive64 + raster500 | 当前默认 | `e891594f` | 已推送 |
| A4 | analytic-gradient、adaptive32 | 淘汰 | `77398635` | 已推送 |
| A5 | raster400 中间值 | 淘汰 | `07fb4c7d` | 已推送 |
| A6 | 已知危险邻域 exact 前置分流 | 接受 | `3e0bef61` | 已推送 |

当前研究分支：`special/approx-safe-research-20260903`。
严格快照分支：`special/exact-fast-validated-20260903`。
严格标签：`bsmpt-v3.1.8-exact-fast-validated-20260903`。

## S0：冻结 exact-fast 严格版本

- 日期：2026-09-03。
- 思路：将此前通过严格 A/B 验证的 CalcGW 实验室状态固化为不可变参照。
- 结果：73+ 严格 A/B 行已验证；严格二进制 SHA256 为
  `b598cc5591bc11aeb378053d4f63e02e18130b3fb865680f008880a05468629e`。
- 结论：作为后续所有非严格试验的 reference，不覆盖、不重建。
- 提交：`f892f02e`，已推送严格分支并创建、推送上述标签。

## A1：analytic-lepton guarded 首遍

- 思路：仅用解析 R2HDM lepton `diff==0` 谱替代对应 Eigen 求解，并用风险检测器
  决定是否回退 exact-fast。
- 样本：NLO 内/外侧、广域 A/B/C、高 SNR Yukawa，共 42 行。
- 结果：首遍接受 6、回退 36；接受行最大 SNR 分量误差 4.36%。发现 A 组状态翻转
  和 B 组 11.58% 总 SNR 反例，均证明不能直接使用无 guard 近似。
- 结论：guard 设计成立，analytic-only 仅保留作消融。
- 提交：`f4ffda0e`，已推送。

## A2：central2 guarded

- 思路：在 analytic-lepton 上把梯度 stencil 换成 central two-point。
- 样本：同一 42 行矩阵，另加 B 假阳性附近四点相对扰动邻域。
- 结果：42 行首遍 2173.759→1597.228 s，减少 26.52%；接受 6 行最大 SNR
  分量误差 3.89%。B 组产生低 SNR 假 GW success；扰动邻域 4/4 双向翻转均被
  状态/非有限值/SNR 地板捕获。
- 结论：只能作为 guarded 首遍；当时升级默认，边缘点不承诺加速。
- 提交：`95ff109e`，已推送。

## A3：adaptive64 与 bounce raster 下探

- 思路：增加 64 点 adaptive exact-solution threshold 搜索，再把 bounce raster
  1000 降到 500；随后用 250、100 探索下限。
- 样本：42 行矩阵、四个 multistep modes、B 假阳性扰动邻域。
- 结果：adaptive64+raster500 首遍 940.724 vs exact-fast 2173.759 s，减少
  56.72%；仍接受 6、回退 36，接受行最大 SNR 分量误差 3.89%。raster250/100
  只再带来约 3–6% 局部收益，NLO 最大误差约 4.78%/5.10%。
- 结论：`central2 + adaptive64 + raster500` 升级为当前 guarded 默认；250/100
  因边际收益小、未覆盖风险扩大而淘汰。
- 提交：`e891594f`，已推送。

## A4：analytic-gradient 与 adaptive32

- 思路：分别用 analytic-gradient 替代 central2，以及把 adaptive grid 64→32。
- 样本：高 SNR type-3/NLO 五点；analytic-gradient 另测 A/B 危险点。
- 结果：analytic-gradient 强信号总 SNR 偏差 0.642%，但慢于 central2；NLO
  五点最大 SNR 分量偏差 13.19%，超过门槛，且首遍 86.047 s 慢于默认
  67.933 s。adaptive32 的非 runtime 输出与 adaptive64 相同，但高 SNR/NLO
  分别为 26.314/96.662 s，均未获益。
- 结论：两项均淘汰，不改默认。
- 提交：`77398635`，已推送。

## A5：raster400 中间值

- 思路：检查 raster500 与 raster250 之间是否有更好折中。
- 样本：NLO 有效边界五点、高 SNR Yukawa 三型、A/B 已知危险点、完整 C 组十点。
- 结果：NLO 最大 SNR 分量偏差 4.77%，66.377 vs 67.933 s；高 SNR 最大偏差
  3.84%，73.367 vs 73.695 s，仅快 0.45%。A/B 风险仍被 guard 捕获。C 组
  10/10 状态/历史一致、最大 SNR 偏差 0.352%，但 238.339 s 慢于默认
  234.463 s。
- 结论：收益不可重复，淘汰，不改默认。
- 提交：`07fb4c7d`，已推送。

## A6：已知危险邻域 exact 前置分流

- 日期：2026-09-03。
- 思路/假设：guard 已知必然回退的点不应先付出近似首遍成本；仅对有保存证据的
  参数锚点/邻域直接运行 exact-fast，不扩大近似接受面。
- 相对基线与唯一变量：保持默认 `central2 + adaptive64 + raster500` 和输出 guard
  不变，只在其前增加 `approx_input_prefilter.py` 与风险锚点表。
- 样本：42 点矩阵全部输入、B 假阳性四点 `-1e-5、-3e-6、+3e-6、+1e-5`
  扰动云；真实集成运行 A 组输入第 10 行。
- 并发数与资源情况：离线判定不运行 CalcGW；集成测试仅一个 CalcGW。
- 状态/history/SNR 检查：42 点中只命中预期 A/B 两点；B 扰动云 4/4 命中；
  A 集成输出除 runtime 外与 exact-fast 参考逐字段一致；高 SNR 安全点不命中。
- exact / 当前默认 / 候选耗时：A 点候选直接 exact 实测 116.945 s；历史默认流程
  会额外运行 41.337 s 近似首遍。B 锚点可避免历史 101.822 s 近似首遍；已验证
  四点扰动云合计可避免 479.302 s 近似首遍。
- 新反例与 guard 是否捕获：无新反例；原 A/B 风险改为在输出 guard 前直接严格算。
- 结论：接受为 safe runner 的保守成本优化。注册表只允许加入已有 exact/approx
  邻域证据的区域，不能外推为通用参数分类器。
- 产物文件：`approx_input_prefilter.py`、`approx_exact_direct_anchors.tsv`、
  `approx_prefilter_a_integration.tsv`、更新后的 safe runner 与报告。
- commit：`3e0bef61`；GitHub：已推送。

## 后续轮次模板

复制以下字段，不得只记录成功项：

```text
## Ax：候选名称
- 日期：
- 思路/假设：
- 相对基线与唯一变量：
- 样本（含输入文件与行号）：
- 并发数与资源情况：
- 状态/history/SNR 检查：
- exact / 当前默认 / 候选耗时：
- 新反例与 guard 是否捕获：
- 结论（升级默认/保留消融/淘汰）：
- 产物文件：
- commit：待提交 / `<hash>`；GitHub：待推送 / 已推送。
```
