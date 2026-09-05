# BSMPT Speed Lab 研究索引

用途：让后续AI按问题找到既有证据，避免从头扫描整个目录或重复失败实验。当前计划仍
以`PROJECT_MASTER_PLAN_ZH.md`和`NON_STRICT_OPTIMIZATION_PLAN_ZH.md`为准。

## 1. 当前基线和入口

| 问题 | 当前事实 | 入口/证据 |
|---|---|---|
| 严格基准是什么 | exact-fast冻结，73+严格A/B行 | `SNAPSHOT.md`、日志S0、`run_calcgw_exact_fast.sh` |
| 当前近似首遍 | central2 + adaptive64 + raster500 + thermal-fast | 日志A2/A3/A8、`run_calcgw_approx_c2_adaptive_r500_thermal_fast.sh` |
| 安全wrapper | 输入风险锚点 + 输出guard + exact fallback | 日志A6、`run_calcgw_approx_safe.sh` |
| 42行性能 | exact 2173.759s，首遍约908.155s，接受6/回退36；回退率不可接受 | `APPROX_SAFE_REPORT_ZH.md`、日志A8 |
| 新的主要缺口 | 状态正常仍可能SNR偏差54.95% | 日志A22、`classified_safe_pilot_counterexample_*` |
| 当前下一步 | N1a定性混淆矩阵与失败阶段账本 | 非严格专项第2、11节；`AGENTS.md`第4节 |

## 2. 已接受并仍在当前路径中的研究

| 轮次 | 内容 | 为什么保留 | 主要产物 |
|---|---|---|---|
| S0 | exact-fast | 严格参照，不得覆盖 | `SNAPSHOT.md`、exact wrapper |
| A2 | central2 | guarded首遍组成；不能裸用 | central2 wrappers、42行结果 |
| A3 | adaptive64+raster500 | 显著降低首遍时间 | adaptive/raster500 wrappers和各组TSV |
| A6 | 已知风险前置exact | 避免必然回退点先付近似成本 | `approx_input_prefilter.py`、anchors TSV |
| A8 | low-temperature fast powers | 当前首遍默认组成 | thermal-fast wrapper、A/B/C/NLO/high-SNR TSV |

所有“接受”仅表示在记录样本和guard条件下保留；不表示完整参数空间已安全。

## 3. 已知危险证据

| 风险 | 事实 | 应查询 |
|---|---|---|
| analytic-lepton状态翻转 | A组bounce结果改变 | 日志A1、`APPROX_SAFE_REPORT_ZH.md`直接近似测试 |
| analytic-lepton SNR超限 | B组总SNR约11.58%偏差 | 同上 |
| central2假GW | B组严格failure而近似success | 日志A2、`central2_false_positive_neighborhood_4*` |
| central2邻域不连续 | 四点扰动均出现双向状态差异 | 非严格报告“假阳性边界邻域” |
| classified E1新反例 | light样本SNR分量偏差54.95% | 日志A22、`classified_safe_pilot_counterexample_*` |
| 历史构建漂移排除 | 当前exact与历史strict仅差0.0281% | `classified_safe_pilot_counterexample_history_vs_current_*` |

N1/N2必须先使上述定性假阳性/假阴性为零；54.95%幅值反例作为双方positive后的次级
数值风险。漏掉任一定性翻转即淘汰候选，不得以最终exact回退掩盖裸近似错误。

## 4. 已淘汰或禁止重复的路线

| 轮次/关键词 | 结论 | 查询位置 |
|---|---|---|
| A4 analytic-gradient | NLO最大SNR偏差13.19%且更慢 | 日志A4、非严格报告 |
| A4 adaptive32 | 输出相同但更慢 | 同上 |
| A5 raster400 | 收益约0.45%且C组更慢 | 日志A5 |
| raster250/100 | 仅增益约3–6%，风险覆盖不足 | 日志A3、非严格报告 |
| A7旧PGO | 画像与当前源码/flags不可信 | 日志A7 |
| A9 high-T fast powers | C组SNR分量偏差43.44% | 日志A9 |
| A10 coefficient pointer | 慢1.41% | 日志A10 |
| A12 thermal cache | 慢0.86% | 日志A12 |
| A13 static log constants | 编译器已折叠 | 日志A13 |
| A14 VEff/mass exact-key cache | 命中率低于门槛 | 日志A14 |
| A15 quark noalias | 慢2.09% | 日志A15 |
| A16 dynamic raster indicator | 无法证明不漏窄峰/零点 | 日志A16 |
| A17 current-source PGO | profile-use缺失画像 | 日志A17 |
| quark 12×12→6×6 | 更慢且改变beta/H | `OPTIMIZATION_REPORT_ZH.md`8.4 |
| 内部多线程旧实现 | 收益小且改变数值 | 严格报告8.5 |
| native vectorization | heap corruption/double-free风险 | 严格报告安全注意 |

重开必须先写“新机制为何不同、以前失败原因为何不再适用、最低停止条件”，不得只换
样本或编译参数重新尝试。

## 5. 数据任务索引

### classified有偏历史数据

- 规则：总纲E1、日志A22/A23；只作为辅助严格参考。
- 抽取：`sample_classified_safe_candidates.py`。
- 当前E1清单：`classified_biased_e1_candidates_240_*`，240点、90 source。
- 比较：`evaluate_classified_approx_pairs.py`。
- 禁止：用`th_passed`筛选安全区、用`th_*`训练router、把空白行当CalcGW失败、从E1
  报告总体成功率。

### 纯BSMPT验证

- E2：围绕safe/unsafe锚点做轴向、随机方向和二分扰动，严格/近似双跑；用于边界。
- E3：规则冻结后按预先定义目标域独立生成；用于接受率、漏放率和端到端性能。
- 当前状态：E2/E3生成器和确认集尚未完成，不得声称已有全空间安全区域。

## 6. 按验证目的找样本

| 目的 | 文件关键词 |
|---|---|
| 42行A/B/C矩阵 | `approx_broad_group_*`、`approx_safe_broad_*` |
| NLO内外边界 | `approx_nlo_*`、`nlo_boundary_*` |
| 高SNR与Yukawa | `approx_high_snr_*`、`stratified_high_snr_*` |
| multistep模式 | `multistepmode*` |
| central2假阳性邻域 | `central2_false_positive_neighborhood_4*` |
| classified E1先导 | `classified_safe_pilot_24_lf_*` |
| 54.95%反例复核 | `classified_safe_pilot_counterexample_*` |
| 新E1辅助池 | `classified_biased_e1_candidates_240_*` |
| exact严格扩展验证 | `full_exact_validation_*`、`extra_validation_*` |

运行前必须从日志确认具体行数、模式和reference，不能仅凭文件名混合不同轮次。

## 7. 工具索引

| 任务 | 工具 |
|---|---|
| 严格运行 | `run_calcgw_exact_fast.sh` |
| 当前guarded运行 | `run_calcgw_approx_safe.sh` |
| 裸近似消融 | `run_calcgw_approx_*`，只能用于配对研究 |
| 同schema逐字段比较 | `compare_outputs.py`、`compare_approx_outputs.py` |
| classified历史配对 | `evaluate_classified_approx_pairs.py` |
| 受控批量执行 | `parallel_calcgw.py`，并发上限2 |
| 单行/少量TSV抽取 | `extract_tsv_rows.py` |
| 判断输出是否应回退 | `approx_needs_fallback.py` |
| 已验证风险输入预筛 | `approx_input_prefilter.py` |

## 8. 研究问题到下一动作

- “为什么回退率高？”：先按guard原因和严格定性类别重建42行账本，不修改算法。
- “如何减少回退？”：可靠positive和可靠fail都应能接受；执行N1/N2，不能只放宽地板。
- “如何判断安全点？”：首先保证正/失败结论一致，再看R0覆盖和C1/C2；幅值为次级。
- “还能否进一步降首遍时间？”：先完成证书；随后按N3顺序逐个做continuation、
  safeguarded shooting、adaptive、SpectrumJet。
- “结果是否可推广？”：只有E3对预先定义的目标分布可推广；E1/E2均不可。
- “下一轮具体做什么？”：读取`AGENTS.md`第4节和非严格专项第11节。
