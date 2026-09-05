# BSMPT speed laboratory

本目录是在不修改原版BSMPT的前提下研究CalcGW加速的隔离实验室。当前项目主目标是
具有严格回退保护的非严格端到端加速；严格快速版的继续优化是次要目标。

## 权威文档

1. `PROJECT_MASTER_PLAN_ZH.md`：唯一项目总纲，定义目标、优先级、边界和证据等级。
2. `NON_STRICT_OPTIMIZATION_PLAN_ZH.md`：主目标的概念架构、具体实现、样本和验收顺序。
3. `RESEARCH_JOURNAL_ZH.md`：逐轮事实、失败路线、产物和Git提交状态。
4. `OPTIMIZATION_REPORT_ZH.md`：已完成严格代码优化的历史技术报告。
5. `APPROX_SAFE_REPORT_ZH.md`：截至A17的非严格实验历史，不再充当当前路线图。
6. `SNAPSHOT.md`：冻结exact-fast版本信息。

其它TSV、JSON、patch和脚本是上述日志引用的实验产物，不应脱离对应轮次解释。

## 不可违反的规则

- 原版BSMPT和主程序不得修改；所有实验只在本目录。
- 最终方案只使用BSMPT/CalcGW，不调用或依赖thdmTools。
- classified是有选择偏差的E1辅助历史数据，不能覆盖全部参数空间，不能定义总体
  安全率或生产路由；任何`th_*`字段禁止成为优化特征。
- 严格路径要求除runtime外逐字段一致。
- 非严格接受结果要求状态、history、有限值结构一致，所有有效SNR分量误差不超过10%。
- unknown、覆盖外、弱信号和证书失败点默认exact-fast。
- 最多两个外层CalcGW；测试内部并行时只运行一个；峰值RSS到18GB停止。

## 当前入口

严格快速版：

```bash
bsmpt_speed_lab/run_calcgw_exact_fast.sh R2HDM input.tsv output.tsv 2 2
```

研究用guarded非严格版（当前仍未获全空间生产认证）：

```bash
bsmpt_speed_lab/run_calcgw_approx_safe.sh R2HDM input.tsv output.tsv 2 2
```

若关心SNR=10和100边界：

```bash
BSMPT_APPROX_SNR_CUTS=10,100 \
  bsmpt_speed_lab/run_calcgw_approx_safe.sh R2HDM input.tsv output.tsv 2 2
```

`run_calcgw_approx_safe.sh`当前首遍是central2 + adaptive64 + raster500 + thermal-fast，
风险点回退到冻结exact-fast。它适合研究，不代表已经证明可用于完整参数空间。

## 当前基线

- 冻结严格分支：`special/exact-fast-validated-20260903`。
- 非严格研究分支：`special/approx-safe-research-20260903`。
- exact-fast已有73+行严格A/B验证。
- 旧42行矩阵：exact-fast 2173.759秒，当前近似首遍约908.155秒；但guard只接受6行、
  回退36行，因此主问题是安全接受率和端到端耗时，而不是继续压缩首遍几个百分点。
- E1辅助先导已确认一个54.95% SNR反例；当前严格复算与历史严格只差0.0281%。

## 当前执行顺序

1. N1：对现有近似增加不改变结果的内部诊断。
2. 用全部已知反例验证诊断的零漏放潜力。
3. N2：action midpoint、raster局部加密、bounce独立复核，分别建立低成本证书。
4. E2：不经过thdmTools的纯BSMPT边界和扰动验证。
5. N3：有证书的continuation、shooting、adaptive和SpectrumJet逐项优化。
6. N4：冻结覆盖/OOD路由，只减少明确不适合近似的浪费。
7. E3：冻结规则后的独立纯BSMPT确认和端到端统计。

严格S1重复gradient消除可穿插进行，但不得挤占主目标长批次资源。

## 数据与比较工具

- `sample_classified_safe_candidates.py`：只凭CalcGW实际执行证据抽取有偏E1辅助集；
  不检查`th_passed`，输出明确标注auxiliary-only。
- `evaluate_classified_approx_pairs.py`：比较历史严格行与新近似结果。
- `compare_approx_outputs.py`：比较同schema严格/近似TSV。
- `parallel_calcgw.py`：受控批量执行；并发不得超过项目资源限制。

E1、E2、E3结果必须分开报告。历史样本通过不能扩大安全区；历史样本失败可以淘汰
候选。只有规则冻结后的E3能够报告明确定义目标分布上的错误率。

## 记录要求

每轮必须记录：假设、唯一变量、数据等级、样本生成方式、并发/RSS、状态/history/SNR、
首遍和端到端耗时、回退率、新反例、结论、文件、commit和push状态。失败实验不得删除
其结论；只删除重复或已被权威文档完整取代的计划说明。
