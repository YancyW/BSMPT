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
| A7 | 复用既有 PGO 构建审计 | 淘汰 | `d5f9afd0` | 已推送 |
| A8 | thermal fast powers | 接受并升级默认 | `7320d3db` | 已推送 |
| A9 | high-temperature fast powers | 淘汰 | `c8edbcaa` | 已推送 |
| A10 | boson coefficient data pointer | 淘汰 | `e016ebaa` | 已推送 |
| A11 | low-temperature fermion/boson split | 保留消融 | `dcd9e6b1` | 已推送 |
| A12 | exact-key thermal last-value cache | 淘汰 | `a9a10e63` | 已推送 |
| A13 | low-temperature static log constants | 淘汰 | `28cfdcc5` | 已推送（状态 `ca15237d`） |
| A14 | VEff/mass exact-key repeat audit | 淘汰 | `396d2bac` | 已推送（状态 `7e68d3f9`） |
| A15 | quark dynamic product noalias | 淘汰 | `e922c704` | 已推送（状态 `bfe13f43`） |
| A16 | dynamic raster safety indicator audit | 停止 | `edc40dbe` | 已推送（状态 `07ba756b`） |
| A17 | current-source mixed-workload PGO | 淘汰 | `af784a96` | 已推送（状态回填待提交） |
| A18 | 选择性PGO只读审计 | 被算法路线取代 | 未实施 | 无独立提交 |
| A21 | 算法级路线重审 | 形成新候选 | `fc9904ca` | 已推送 |
| A22 | classified 24点先导 | 口径部分作废，反例有效 | `d30145f1` | 已推送 |
| A23 | 有偏数据口径纠正 | 接受 | `f1172400` | 已推送 |
| A24 | 文档体系重构 | 当前权威结构 | `ee7869cb` | 已推送 |

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

## A7：复用既有 PGO 构建审计

- 日期：2026-09-03。
- 思路/假设：在不改变数值算法的情况下，将现有 PGO binary 叠加到 guarded 首遍，
  争取历史测试中约 2–5% 的编译级收益。
- 相对基线与唯一变量：计划仅替换 `BSMPT_CALCGW_BINARY`，其余保持 A6 默认。
- 样本：计划使用高 SNR 与 NLO 边界；构建资格审计未通过，因此未启动 CalcGW。
- 并发数与资源情况：零个 CalcGW。
- 审计结果：`build-pgo-gcc` 的 cache 仍为 `-fprofile-generate`，binary 是 158 MB
  instrumentation 版本而非 profile-use 成品；目录有 35 个 `.gcda`，但源码时间戳
  晚于训练 binary，不能证明画像与当前源码匹配。构建目录为 542 MB。
- 结论：禁止直接用于测速或 safe runner。若未来重启 PGO，必须从当前源码重新生成
  覆盖广域/NLO/强信号的训练集，并在隔离的新目录构建 profile-use 版本。
- 产物文件：仅本审计日志，无新 binary、无 CalcGW 输出。
- commit：`d5f9afd0`；GitHub：已推送。

## A8：thermal fast powers

- 日期：2026-09-03。
- 思路/假设：低温 thermal 展开中大量整数/半整数 `pow` 可用乘法与一次 `sqrt`
  代替，在非严格模式允许的舍入差异内减少约 200 万次 VEff 调用的成本。
- 相对基线与唯一变量：`central2 + adaptive64 + raster500` 不变；独立
  `build-approx-thermal` 仅增加 `BSMPT_USE_THERMAL_FAST_POWERS=1`。严格 binary
  不重建，SHA256 仍为 `b598cc...629e`。
- 样本：42 行 NLO/A/B/C/高 SNR 矩阵，两次 type-3 同负载配对，四种 multistep
  modes，A/B 已知危险点，端到端 safe 高信号点。
- 并发数与资源情况：构建 `-j2`；CalcGW 最多两个并发。
- 状态/history/SNR 检查：接受面仍为 A 三行和高 SNR 三行。A 接受行相对严格
  最大 SNR 分量误差 0.026%，高 SNR 最大 3.74%；四种 mode 最大 3.80%。C 组
  10/10、NLO 有效 5/5 状态/history 一致，NLO 外侧 4/4 保持拒绝。
- exact / 当前默认 / 候选耗时：42 行 2173.759 / 940.724 / 约 908.155 s；候选
  相对旧默认再降 3.46%，相对 exact 降约 58.22%。C 组 234.463→220.950 s，
  高 SNR 三型 73.695→71.443 s。
- 新反例与 guard 是否捕获：没有新反例。A/B 只复现已知状态差异；A6 会在首遍前
  直接 exact。B 锚点候选 SNR 对旧默认可变 11.21%，但不会进入候选或被接受。
- 结论：接受并升级 guarded 默认首遍；严格 fallback 继续使用未修改 exact binary。
- 产物文件：thermal 源码开关、独立 wrapper、control/candidate/mode/integration TSV。
- commit：`7320d3db`；GitHub：已推送。

## A9：high-temperature fast powers

- 日期：2026-09-03。
- 思路/假设：在 A8 之上复用 `sqrt(x)`、`1/sqrt(x)` 与 `1/x`，替代高温展开
  `JInterpolatedHigh` 中 `pow(x, ±l/2)`，进一步减少 thermal 成本。
- 相对基线与唯一变量：A8 low-temperature fast powers 保持开启，只增加默认关闭的
  `BSMPT_USE_THERMAL_FAST_HIGH_POWERS=1`。
- 样本：高 SNR type-3 同负载配对、NLO 有效五点同负载配对、完整 C 组十点。
- 并发数与资源情况：增量构建 `-j2`；CalcGW 最多两个并发。
- 状态/history/SNR 检查：三组状态/history 均不变；高 SNR 相对 A8 新增最大 SNR
  偏差 0.208%，NLO 相对严格最大 5.07%；但 C 组第 9 行 SNR 分量相对严格偏差
  43.44%，超过 10% 硬门槛。
- exact / A8 / 候选耗时：高 SNR type-3 23.133→22.338 s；NLO 五点
  66.585→66.782 s，无收益；C 组 220.950→209.593 s，快 5.14%。
- 新反例与 guard 是否捕获：C 第 9 行是新的数值敏感反例；当前因弱信号会回退，
  但其大幅偏差说明不能向未覆盖接受区推广。
- 结论：淘汰，不加入 A8 wrapper；实验开关默认关闭并保留反例供反查。
- 产物文件：thermal high source switch、high/NLO/C control/candidate TSV。
- commit：`c8edbcaa`；GitHub：已推送。

## A10：boson coefficient data pointer

- 日期：2026-09-03。
- 思路/假设：A8 boson 低温快路径每次读取二、三阶系数，尝试一次取得连续 data
  pointer，替代两次 `GetCoefficentAtOrder`。
- 相对基线与唯一变量：A8 不变，仅增加默认关闭的
  `BSMPT_USE_THERMAL_COEFF_DATA=1`。
- 样本：高 SNR type-3 同负载配对。
- 并发数与资源情况：增量构建 `-j2`；两个 CalcGW 并发。
- 状态/history/SNR 检查：所有非 runtime 字段逐位一致。
- exact / A8 / 候选耗时：A8 control 23.317 s，候选 23.646 s，慢 1.41%。
- 新反例与 guard 是否捕获：无数值反例；纯性能失败。
- 结论：淘汰，不加入 A8 wrapper；原 accessor 已内联，额外环境分支无收益。
- 产物文件：source experiment、`thermal_coeff_control_high.tsv`、
  `thermal_coeff_candidate_high.tsv`。
- commit：`e016ebaa`；GitHub：已推送。

## A11：low-temperature fermion/boson split

- 日期：2026-09-03。
- 思路/假设：拆分 A8 的 fermion 与 boson 快幂，定位速度和数值偏差来源，并寻找
  更稳健的子集。
- 相对基线与唯一变量：新增默认关闭的 `BSMPT_USE_FERMION_FAST_POWERS` 与
  `BSMPT_USE_BOSON_FAST_POWERS`；原 A8 组合开关语义不变。
- 样本：高 SNR type-3 初筛；fermion-only 与无 fast 再做同负载复测。
- 并发数与资源情况：增量构建 `-j2`；两个 CalcGW 并发。
- 状态/history/SNR 检查：fermion-only 相对无 fast 所有非 runtime 字段逐位一致；
  boson-only 改变总 SNR 与 beta/H。A8 的主要数值偏差来自 boson 部分。
- exact / control / 候选耗时：复测 control 24.022 s，fermion-only 23.711 s，快
  1.29%；首轮 fermion/boson 分别 23.161/23.018 s，不跨轮直接比较。
- 新反例与 guard 是否捕获：无新状态反例；boson-only 的数值变化说明不能将其当作
  严格子优化。
- 结论：fermion-only 保留为更保守消融，不替换完整矩阵已验证且收益更高的 A8；
  boson-only 不单独推广。
- 产物文件：split source switches 与四个 type-3 TSV。
- commit：`dcd9e6b1`；GitHub：已推送。

## A12：exact-key thermal last-value cache

- 日期：2026-09-04。
- 思路/假设：由 Luna 只读审计设计默认关闭的连续 bitwise-identical thermal 调用
  observer；命中率超过 5% 后，实现 thread-local last-value exact-key cache。
- 相对基线与唯一变量：A8 不变；observer 用
  `BSMPT_PROFILE_THERMAL_REPEATS=1`，缓存用默认关闭的
  `BSMPT_USE_THERMAL_EXACT_LAST_CACHE=1`。键包含 kind、diff 与 x 原始位模式，无
  容差/量化，缓存不跨线程。
- 样本：高 SNR type-3 profiling 和同负载 control/candidate。
- 并发数与资源情况：增量构建 `-j2`；profiling 一个 CalcGW，配对两个 CalcGW。
- 命中率：boson d0 `1,345,263/23,133,504`（5.82%），fermion d0
  `12,804,236/40,483,632`（31.63%），合计约 22.3%。
- 状态/history/SNR 检查：cache control/candidate 所有非 runtime 字段逐位一致。
- exact / control / 候选耗时：control 21.970 s，cache 22.158 s，慢 0.86%。
- 新反例与 guard 是否捕获：无数值反例；TLS、位键与分支成本超过被省 thermal
  算术，属于纯性能失败。
- 结论：淘汰，不扩展完整矩阵；observer 与默认关闭的实验开关保留供反查。
- 产物文件：profiler/thermal instrumentation、profile/control/candidate TSV。
- commit：`a9a10e63`；GitHub：已推送。

## A13：low-temperature static log constants

- 日期：2026-09-04。
- 思路/假设：由 Luna 只读审计提出，将 fermion/boson 低温展开中只依赖常量的
  `cf/cb` 改成函数内 `static const`，期望避免重复常数 `log`。
- 相对基线与唯一变量：只做当前 A8 `build-approx-thermal` 目标文件的反汇编审计；
  没有修改任何源码、wrapper、binary 或默认路径。
- 样本：不进入 CalcGW 数值阶梯；最低阶梯已能判定候选不存在实际工作量。
- 并发数与资源情况：未运行 CalcGW，未构建；严格 binary SHA256 保持不变。
- 状态/history/SNR 检查：无代码或 binary 变化，因此不产生新的数值结果。
- 性能证据：当前编译参数含 `-O3`；目标文件中 `cf/cb` 已由 `.LC` 只读常量载入，
  `JfermionInterpolatedLow4Exact`/`JbosonInterpolatedLow` 热路径的 `log` relocation
  只对应运行时变量 `x`。手工静态化不会消除一次现存的常数 `log`。
- 新反例与 guard 是否捕获：不适用。
- 结论：最低阶梯淘汰；编译器已完成该优化，不实施无效源码改写，也不消耗配对运行。
- 产物文件：本日志与路线图状态；无新 TSV。
- commit：`28cfdcc5`；GitHub：已推送（状态回填 `ca15237d`）。

## A14：VEff/mass exact-key repeat audit

- 日期：2026-09-04。
- 思路/假设：由 Luna 只读审计后，为 VEff、HiggsMasses、QuarkMasses 增加默认
  关闭的紧邻 exact-bit key observer；先测命中率，再决定是否实现 last-value cache。
- 相对基线与唯一变量：A8 binary 数学路径不变；仅 profile 进程设置
  `BSMPT_PROFILE_VEFF_REPEATS=1`。key 含 model 地址、全部场值原始位、温度、diff
  和 order，三类分别维护 thread-local 上一条记录；禁止容差或量化。
- 样本：高 SNR type-3 同负载 control/profile；命中率在最低运行阶梯已低于停止线，
  因而按路线图不扩展 NLO/C 长批次。
- 并发数与资源情况：增量构建 `-j2`；两个 CalcGW 并发，无其他 CalcGW。
- 状态/history/SNR 检查：control/profile 一行除 runtime 外逐字段零容差一致。
- 命中率：VEff `789/2,006,315`（0.039%），Higgs `1/53`（1.89%），Quark
  `801/2,006,317`（0.040%），全部低于 3% 停止线。
- 性能证据：control 27.315 s，profile 27.650 s；observer 慢 1.23%。此数值只用于
  确认 instrumentation 成本，不将并发单次计时当成候选性能结论。
- 清理复核：移除 observer 并重建后的同一点为 26.561 s，和 profile 前 control
  除 runtime 外逐字段一致；单次绝对耗时仍受同机负载影响，不跨轮宣称 0.754 s 收益。
- 新反例与 guard 是否捕获：无数值反例；本轮没有实现缓存。
- 结论：淘汰 VEff/Higgs/Quark 紧邻 exact-key last-value cache；TSV 与命中率结论
  供反查，不扩展更复杂窗口缓存。observer 在测量后从默认构建源码移除，避免每次
  VEff/质量调用的关闭判断污染正式首遍性能。
- 产物文件：提交历史中的临时 profiler/model instrumentation、
  `p3_high_snr_type3_input.tsv`、
  `p3_repeat_control_high.tsv`、`p3_repeat_profile_high.tsv`。
- commit：`396d2bac`；GitHub：已推送（状态回填 `7e68d3f9`；observer 清理
  `e922c704`）。

## A15：quark dynamic product noalias

- 日期：2026-09-04。
- 思路/假设：由 Luna 审计 P4 后，对 quark 动态路径的
  `MassMatrix = MIJ.conjugate() * MIJ` 单独测试 Eigen `noalias()` 赋值，尝试减少
  临时对象或别名检查，不改变矩阵维度、乘法表达式、solver 或本征值处理。
- 相对基线与唯一变量：A8 及所有既有开关不变；候选仅设置默认关闭的临时
  `BSMPT_USE_QUARK_MASS_PRODUCT_NOALIAS=1`。fixed12 的 diff<=0 路径不受影响，
  动态导数路径可覆盖该候选。
- 样本：高 SNR type-3 同负载 control/candidate；最低性能阶梯失败后停止。
- 并发数与资源情况：增量构建 `-j2`；两个 CalcGW 并发。
- 状态/history/SNR 检查：一行除 runtime 外逐字段零容差一致。
- exact / control / 候选耗时：control 28.256 s，candidate 28.846 s，候选慢
  2.09%；未达到 1% 收益门槛。
- 新反例与 guard 是否捕获：无数值反例；纯性能失败。
- 结论：淘汰并从当前源码删除实验开关；历史 6x6、fixed eigensolver、stack
  matrix、pair reserve 等均已失败或不适用，P4 矩阵临时/容量路线判定耗尽。
- 产物文件：`p4_noalias_control_high.tsv`、`p4_noalias_candidate_high.tsv`；实验源码
  只留在本轮工作历史，不进入默认 binary。
- commit：`e922c704`；GitHub：已推送（状态回填 `bfe13f43`）。

## A16：dynamic raster safety indicator audit

- 日期：2026-09-04。
- 思路/假设：由 Luna 只读审计 raster500/1000 路径，寻找能在首遍运行时识别必须
  升级 raster1000 的内部指标，同时要求对 A/B/C/NLO 边缘零漏检且不扩大接受面。
- 相对基线与唯一变量：仅审计 `action_calculation.cpp` 中 rasterized dV/dl、精确
  `Calc_dVdl` 和 RK5 消费路径；没有修改源码、wrapper 或 binary。
- 样本：复核既有 raster 降阶结论和验证要求；因不存在可证明安全的指标，不启动
  新 CalcGW 批次。
- 并发数与资源情况：零个 CalcGW，未构建；严格 binary 不变。
- 状态/history/SNR 检查：没有运行路径变化，不产生新数值结果。
- 审计结果：可诊断的候选是在部分 RK5/转折/端点位置比较 raster500 插值与精确
  dV/dl 的归一化残差；但有限采样会漏掉未采样窄峰、零点移动、RK5 未经过区间及
  后续 SNR 放大。阈值取零才可能零漏检，却等价于所有点回到 raster1000。
- 新反例与 guard 是否捕获：现有有限 paired 数据只能校准已知点，不能证明高维参数
  空间无漏检；把诊断阈值用于接受会扩大未经证明的接受面。
- 结论：P5 停止，不实施动态 raster；保留固定 raster500 + 输出 guard + exact
  fallback 的现状。
- 产物文件：本日志与路线图状态；无新 TSV。
- commit：`edc40dbe`；GitHub：已推送（状态回填 `07ba756b`）。

## A17：current-source mixed-workload PGO

- 日期：2026-09-04。
- 思路/假设：旧 PGO 画像与当前源码/flags 不匹配，故从当前 overlay 新建隔离
  control/generate 树，用 A/B/C、NLO 内外和高 SNR 混合训练后再构建 profile-use。
- 相对基线与唯一变量：control 与 generate 使用当前 A8 build 的 `/usr/bin/gcc/g++`
  14.3.1、Release、相同 `-march=nocona -mtune=haswell -O2/-O3` 等 flags；generate
  只增加定向 `-fprofile-generate`。严格树与严格 binary 均未重建。
- 样本：`pgo_current_training_mixed_6.tsv` 六行，覆盖 broad A/B/C、NLO valid/
  invalid 和高 SNR type-3；单个插桩 CalcGW 串行训练，输出 6/6 行。
- 并发数与资源情况：control/generate 构建各 `-j1` 同时运行，总编译并发 2；训练
  仅一个 CalcGW；profile-use 构建 `-j2`，首个失配警告后立即终止。
- 画像结果：训练前 `.gcda=0`，训练后生成 35 个当前路径画像。旧 Conda GCC 14.3.0
  配置因当前 glibc 版本符号无法链接而弃用，未复用旧 35 个画像。
- 资格失败：在同一 generate 目录切换 `-fprofile-use -fprofile-correction` 后，首轮
  构建报告 `ClassPotentialOrigin_deprecated.cpp.gcda profile count data file not found`。
  按路线图“任一缺失/checksum/profile mismatch 即停止”的规则，不生成或测试部分
  profile-use binary。
- 状态/history/SNR 检查：训练输出仅用于画像；profile-use 未通过构建资格，故不进入
  数值/性能验收，也不改变 guarded 默认。
- 新反例与 guard 是否捕获：无数值反例；构建画像覆盖失败。
- 结论：淘汰本轮 PGO，不绕过 missing-profile 警告，不复用旧画像。
- 产物文件：`pgo_current_training_mixed_6.tsv`、训练输出；build 与 `.gcda` 由
  `.gitignore` 排除。
- commit：`af784a96`；GitHub：已推送，状态回填提交随后推送。

## A18：选择性热点翻译单元 PGO（未实施，已被取代）

- 日期：2026-09-04。
- 思路/假设：A17 的全程序 profile-use 因未执行翻译单元缺少 `.gcda` 而失败；只对
  已被混合 R2HDM 训练实际覆盖、且属于 CalcGW 热路径的翻译单元施加 generate/use，
  其余对象保持与 A8 control 相同 flags，可避免部分画像 binary 和 missing-profile。
- 相对基线与唯一变量：保持 A8 数学路径、central2、adaptive64、raster500、编译器
  和依赖不变，只研究逐源文件 PGO flags；严格 binary 不重建。
- 固定资格条件：先核对35个 `.gcda` 与目标对象映射；任何已选择热点对象缺失、
  checksum mismatch 或 profile warning 即淘汰。通过后才进入高 SNR 配对测试。
- 并发限制：构建最多 `-j2`，CalcGW 最多两个。
- 结论：只读审计后未实施；项目重点已经转向非严格多保真数值证书，不再列入当前队列。
- commit：无独立实现提交。

### A18 Luna 并行审计结论

- 三个只读任务分别审计 PGO机制、热点集合及 A19/A20 后备路线；均未编辑、构建或
  运行 CalcGW。
- 35个 current-system `.gcda` 来自 GCC 14.3.1。deprecated Origin函数当前无调用，
  因而没有画像；四个最小化器源文件又被多个 target 重复编译，不能全局启用PGO。
- M0固定为8个有画像、无重复对象歧义的热点翻译单元，约覆盖累计热点排序的88%；
  预期收益1%–4%，仍执行2%硬门槛。
- 实现应使用按源文件、`-o`对象路径和target三重匹配的实验室 compiler-launcher；
  新建单一 generate/use 构建树，并把 missing-profile、coverage-mismatch 升为错误。
- A19后备优先级为 lld safe ICF、干净 non-PIE，再考虑单文件编译选项；A20优先批量
  guarded runner，保持逐行 guard 与 exact fallback，预期在大量廉价点时收益最大。

## A21：算法级路线重审

- 日期：2026-09-05。
- 思路：停止把个位数编译优化作为主路线，从真实调用次数、数据依赖和残差证书重新
  设计严格与guarded算法。
- 分工：三个Luna只读审计分别覆盖bounce温度/profile continuation、R2HDM谱与导数、
  确定性并行与phase tracing；主代理核对现有源码和历史，统一排序。
- 核心账本：严格代表点bounce约26.9秒，其中24次Solve1DBounce约21.5秒；每action
  的1001点raster约16016次VEff，dense threshold约41041次VEff。
- 新结论：路径几何warm start已经存在；真正未复用的是径向profile、射击括区间和
  最佳l0。最大严格机会是重复求值消除、threshold/raster固定索引并行；最大guarded
  机会是SpectrumJet、profile continuation、误差驱动adaptive采样。
- 产物：原路线图内容已整理进 `PROJECT_MASTER_PLAN_ZH.md` 和
  `NON_STRICT_OPTIMIZATION_PLAN_ZH.md`；旧文件可从Git历史恢复。
- 下一实验：S1 `LocateMinimum` 同点重复gradient消除；待实现和验证。
- commit：`fc9904ca`；GitHub：已推送。

## 路线图完成后的残余审计

- 日期：2026-09-04；由 Luna 对 A1–A17、失败路线与当前 overlay 做只读复核。
- 结论：当前局部优化空间耗尽，未发现仍满足“单变量、预期总收益至少 2%、边缘
  安全”的新候选。thermal/VEff cache、常量折叠、矩阵临时、动态 raster、旧 PGO、
  LTO、固定降 raster、analytic gradient 等均已有直接失败证据或违反安全前提。
- 只有两类新信息值得重开：一是能覆盖全部链接对象、零 missing/mismatch 警告的
  当前源码 PGO 画像；二是新编译器/Eigen 实现并附算法等价证明。二者仍须从固定
  A/B/C/NLO/高 SNR 阶梯重新验收，不能沿用旧结论。

## A22：classified 条件安全区域与 24 点先导

- 日期：2026-09-05。
- 思路/假设：只在 thdmTools 已通过的条件样本中学习近似安全性；`neither` 不作
  CalcGW 负例，使用历史严格输出降低验证成本。
- 样本：三台机器 `both`/`thdmtool_only` 只读流式筛选；生成 240 点、100 source
  的候选清单，并以其中24点作接口和配对先导。
- 并发数与资源情况：始终一个 CalcGW；未扫描绝大多数 `neither`，未写 MultiNest。
- 接口纠错：首遍 CRLF 破坏 CalcGW 拼接输出，已强制 LF 并同样本重跑；错误输出
  不进入统计。
- 状态/history/SNR 检查：LF 先导23/24满足配对条件；19个弱/零信号只可严格路由；
  4个有信息量接受点最大误差分别约0.015%、0.45%、5.29%、0.043%。
- 新反例：light物理边界一点偏差54.95%。当前 exact-fast 与历史严格仅差0.0281%，
  确认是近似反例而非版本漂移。
- 结论：现阶段只能确定候选域，不能上线区域路由；light边界、弱信号、多步/异常、
  Yukawa II--IV保持严格。下一步扩大配对并按source/机器留出验证。
- 产物：旧计划已删除、事实保留于本轮日志；`sample_classified_safe_candidates.py`、
  `evaluate_classified_approx_pairs.py`、`classified_safe_pilot_24_*`、
  `classified_safe_candidates_240_lf_*`。
- commit：`d30145f1`；GitHub：已推送。

## A23：纯BSMPT规约与有偏数据口径纠正

- 日期：2026-09-05。
- 思路/假设：classified仅作为经过上游选择的E1辅助历史集；优化、router和证书只能
  使用BSMPT输入及本次BSMPT诊断，安全边界必须由无thdmTools的E2/E3确认。
- 纠错：A22中以`th_passed=true`筛选及把历史密集分支解释为候选安全域的口径作废；
  24点数值配对和已由当前exact复核的54.95%反例仍是有效E1证据。
- 工具：`sample_classified_safe_candidates.py`现只凭完整BSMPT输入、runtime/核心状态和
  source认证CalcGW实际执行；默认扫描both/gw_only/thdmtool_only，neither因体量巨大
  仅显式请求时扫描，且类别只用于文件定位，不能作为特征。
- 新E1清单：240点，90个source，office/qiushi/yuanyuan为83/80/77；resolved signal
  35、SNR边界75、数值边界75、严格弱/无信号55。该分布不得用于总体比例推断。
- 优化设计：严格优先S1重复gradient、S2/S3固定索引并行；guarded优先G1 coarse/fine
  bounce后验证书，再做continuation、safeguarded射击和adaptive节点。
- 统计：E1只能淘汰和找反例；E2决定候选局部边界；只有预先冻结、无筛选的E3可报告
  明确定义目标分布上的漏放率及端到端收益。
- 产物：规约已整理进 `PROJECT_MASTER_PLAN_ZH.md` 和
  `NON_STRICT_OPTIMIZATION_PLAN_ZH.md`、修正后的抽样器、
  `classified_biased_e1_candidates_240_*`。
- commit：`f1172400`；GitHub：已推送，状态回填提交随后推送。

## A24：文档体系重构

- 日期：2026-09-05。
- 原因：README、三个路线图、A17后续计划、综合规约及错误classified计划职责重叠，
  导致主目标和当前执行顺序不清楚。
- 新结构：`PROJECT_MASTER_PLAN_ZH.md`为唯一总纲；
  `NON_STRICT_OPTIMIZATION_PLAN_ZH.md`为主目标唯一执行方案；README仅导航。
- 主/次目标：非严格端到端加速为主，重点降低安全识别成本和过度回退；严格S1--S3
  为次要目标，不占用主目标长批次资源。
- 保留：逐轮日志、严格优化报告、非严格实验报告、严格快照，确保事实可反查。
- 删除：已被吸收的多个路线图、后续计划、错误口径计划，以及已经完成/淘汰的孤立
  设计提案。内容仍可从Git历史恢复，相关实验事实仍在报告和日志中。
- 下一执行：N1只读内部诊断，不改变近似数值路径；先验证54.95%等已知反例能否被
  低成本识别，再决定N2双分辨率证书。
- commit：`ee7869cb`；GitHub：已推送，状态回填提交随后推送。

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
