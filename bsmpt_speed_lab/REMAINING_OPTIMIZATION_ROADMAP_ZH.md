# CalcGW guarded 剩余优化路线图

更新日期：2026-09-04。本文件在实验前固定方向、优先级、验证阶梯与停止条件；
逐轮事实和提交状态继续写入 `RESEARCH_JOURNAL_ZH.md`。

## 当前基线

- 严格参照：分支 `special/exact-fast-validated-20260903`、提交 `f892f02e`，禁止修改。
- guarded 分支：`special/approx-safe-research-20260903`。
- 当前首遍：central2 + adaptive64 + raster500 + A8 low-temperature fast powers。
- A6 对已验证 A/B 风险锚点和 B 扰动邻域直接 exact-fast。
- 42 行：exact-fast 2173.759 s；旧 raster500 940.724 s；A8 约 908.155 s。
  guard 接受 6 行、回退 36 行，接受行已测最大 SNR 分量偏差 3.80%。

## 不变约束

- 所有编辑、构建和输出只允许在 `bsmpt_speed_lab/`。
- 严格 binary SHA256 必须保持
  `b598cc5591bc11aeb378053d4f63e02e18130b3fb865680f008880a05468629e`。
- 最多两个 CalcGW，编译最多 `-j2`；长批次前检查现有进程。
- 子 agent 只读审计，不运行 CalcGW、不构建、不编辑、不做 Git 操作。
- 主 agent 统一实施和验证，每次只改变一个独立变量。
- 成功和失败轮次都保存关键 TSV、更新日志、提交并推送。

## 验收门槛

1. guard 接受行的全部 status 和 `transition_history` 与严格参照一致。
2. 接受行每个 SNR 分量相对严格参照不超过 10%；7% 以上进入警戒复验。
3. A/B 边缘和 B 扰动云必须被前置 exact 或输出 guard 捕获。
4. NLO 外侧保持拒绝；NLO 内侧、C 组和四种 multistep mode 不新增状态翻转。
5. 单点收益至少重复两次，或不少于五点的组总耗时至少改善 2%。
6. 任一边缘/弱信号 SNR 分量偏差超过 10% 默认淘汰，即使当前会回退；只有严格、
   不扩大接受面的前置 exact 规则能完整隔离时才可重新评估。
7. exact fallback 始终使用严格 binary，近似环境变量不得泄漏。

## 固定验证阶梯

候选在最低阶梯失败即停止：

1. 历史重复检查、静态审计、编译、语法、严格 binary 哈希。
2. 高 SNR type-3 同负载 control/candidate 两次。
3. NLO 有效五点和 A/B 两个危险点。
4. 完整 C 十点和高 SNR Yukawa 2/3/4。
5. 完整 A/B，并按实际 guard 分开统计接受/回退行。
6. multistepmode 0/1/2/auto 和 B 四点扰动云。
7. 42 行汇总及端到端 safe 接受/前置 exact 集成。

## 候选队列

### P1：exact-key thermal result cache

状态：A12 已完成并淘汰。type-3 连续 exact-key 命中约 22.3%，但缓存慢 0.86%。

先用纯计数 instrumentation 测量连续或小窗口内 `(kind,diff,x-bit-pattern)` 完全相同
的调用命中率。命中足够高才实现 thread-local 小缓存；禁止量化、容差和近邻复用。
潜在收益中等，数值风险低。命中率低于 5%，或两次配对收益低于 1% 即停止。

### P2：A8 低温公共子表达式

状态：A13 已完成静态审计并淘汰。`cf/cb` 的常数对数表达式已被当前 `-O3`
编译器折叠为只读常量，目标文件热路径中的 `log` 调用只剩变量 `x`；无需再做等价改写。

只测试不改变展开阶数的局部算术，每项独立开关。A10 coefficient pointer 已无收益，
不得重复。潜在收益低、舍入风险中；边缘 SNR 超 10% 或五点收益低于 2% 即停止。

### P3：VEff/质量谱 exact-key 命中率审计

状态：A14 已完成并淘汰。高 SNR type-3 的紧邻 exact-bit 命中率为 VEff 0.039%、
Higgs 1.89%、Quark 0.040%，均低于 3% 停止线；不实现结果缓存。

先统计同一 model instance 内场值、温度和 diff 逐位相同的重复调用，再决定是否为
quark/Higgs 增加 last-value 小缓存。必须覆盖所有状态，生命周期不得跨参数点。
quark/Higgs 热点约 8.68/6.57 s，潜在收益高；命中率低于 3% 或键无法证明即停止。

### P4：矩阵构造公共子表达式与容量

只优化填充和临时内存，不改变 Eigen 求解器、矩阵维度、本征值复制或求和顺序。
历史 quark 6x6、fixed eigensolver、stack matrix 已失败。配对收益低于 1% 即停止。

### P5：动态 raster（延后）

只有找到能在 A/B/C/NLO 边缘零漏检 raster500/1000 差异的内部指标才实施；不得
仅用 bounce 成败判断。指标失败立即停止。

### P6：全新 PGO（最后）

旧画像失配。算法候选耗尽后，才从当前源码重新 generate，用 A/B/C/NLO/高 SNR
混合训练并在隔离目录 profile-use。画像警告、flags 不匹配或收益低于 2%即停止。

## 禁止重复路线

raster400/250/100 固定下调、adaptive32、analytic-gradient、high-temperature fast
powers、boson coefficient pointer、quark 12x12→6x6、directional dV/dl、RK5
fixed/workspace、旧 PGO 复用、容差键缓存、未经验证扩大 A6 风险半径，均不再测试。

## Luna 分工与用量控制

1. 同时只运行一个 Luna；每次只审计一个候选。
2. Luna 只交付源码落点、历史重复检查、风险和最小实验，不做实际运行。
3. 主 agent 每次只选择一个候选实施，并统一控制最多两个 CalcGW。
4. 最低阶梯失败立即终止该候选，避免浪费完整矩阵用量。
5. 每轮立即更新日志中的思路、结果、结论、commit 和 GitHub 状态。
