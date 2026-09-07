# Region-v5 type-I专项E3协议

预注册日期：2026-09-07

从本轮起，非严格优化的新增研究只聚焦Yukawa type-I。既有type-II/III/IV结果保留为
历史证据，但不再扩样、调参或占用优化资源。本协议在查看新样本结果前冻结，结果不得
用于修改v5认证半径。

## 固定样本

全部`yuktype=1`，使用与建域阶段不同的随机种子，每个区域2点：

| 组 | 来源 | delta | seed | 路由预期 |
|---|---|---:|---:|---|
| T1 | broad-A row1 positive | `5e-5` | 61007 | 正点认证 |
| T2 | broad-A row2 positive | `5e-5` | 61107 | 正点认证 |
| T3 | broad-C row3 positive | `8e-5` | 61207 | 正点认证 |
| T4 | broad-C row5 positive | `5e-5` | 61307 | 正点认证 |
| T5 | broad-A row4 coexistence | `5e-3` | 61407 | coexistence认证 |
| T6 | broad-B row1 coexistence | `5e-5` | 61507 | coexistence认证 |
| T7 | broad-C row1 coexistence | `5e-5` | 61607 | coexistence认证 |
| T8 | broad-B row10 bounce | `2e-5` | 61707 | 域外严格回退 |
| T9 | broad-A row9 late failure | `4e-5` | 61807 | 域外严格回退 |

共18点。T1/T2/T5/T9来自`stratified_broad_group_a_10.tsv`；T6/T8来自B组；
T3/T4/T7来自C组。统一使用`generate_anchor_neighborhood.py --random-directions 2
--max-points 2`。

## 验收条件

1. strict与v5正/失败双向一致，FP=FN=0。
2. T1--T7必须14/14由预期证书实际接受，回退不计通过。
3. T8/T9必须4/4严格回退。
4. 每个被认证组加入NLO前缀后必须有正向收益；全18点总降时至少20%。
5. 正点报告最大SNR分量误差；失败阶段漂移单独报告。
6. 本批只说明预注册type-I分层，不外推到完整type-I连续参数空间。

## 冻结后结果

- 18点得到TP=8、TN=10、FP=FN=0；T1--T7共14点全部由预期证书接受，T8/T9共4点
  全部严格回退。
- T1--T7修正降时依次为58.50%、61.79%、65.42%、60.13%、30.92%、25.70%、
  27.06%。全批严格`1553.213 s`，v5加NLO前缀`831.412 s`，降低46.47%。
- 可比较SNR字段中的最大分量误差为0.384%。但T3两点严格/近似相变历史schema均不一致，
  因此虽然主契约通过，仍将该域从安全近似中撤销并由region-v6前置严格。
- region-v6复验T3两点均direct exact；与严格输出除runtime外逐字段完全一致。
