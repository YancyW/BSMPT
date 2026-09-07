# Region-v5 独立E3确认协议

预注册日期：2026-09-07

本文件在查看本批严格/候选结果前冻结。E3结果不得用于调整region-v5认证半径、锚点或
决策规则；若认证域出现FP/FN，撤销对应认证域并将该轮记为失败。

## 固定样本

采用与E20--E23不同的随机种子和扰动幅度：

| 组 | 来源 | 半径 | 点数 | seed | 预期路由用途 |
|---|---|---:|---:|---:|---|
| C1 | high-SNR Yukawa type-2 | `8e-5` | 4 | 52007 | v5正点认证域 |
| C2 | high-SNR Yukawa type-4 | `2e-5` | 4 | 52107 | v5正点认证域 |
| C3 | broad-B row1 coexistence | `8e-5` | 4 | 52207 | v5失败认证域 |
| C4 | broad-C row1 coexistence | `5e-5` | 4 | 52307 | v5失败认证域 |
| C5 | broad-B row10 bounce | `3e-5` | 2 | 52407 | 域外严格回退控制 |
| C6 | broad-A row9 late failure | `6e-5` | 2 | 52507 | 域外严格回退控制 |

总计20点。生成器固定为`generate_anchor_neighborhood.py --random-directions N
--max-points N`。C1/C2分别来自`stratified_high_snr_yukawa_3.tsv`第1/3行；C3/C5来自
`stratified_broad_group_b_10.tsv`第1/10行；C4来自`stratified_broad_group_c_10.tsv`
第1行；C6来自`stratified_broad_group_a_10.tsv`第9行。

预注册勘误：在生成任何样本、查看任何结果之前，静态检查发现生成器使用
`delta*max(abs(value),1)`，而router使用`radius*abs(reference)`。因此将C1--C4的
生成delta修正为上表数值，确保小参数方向也位于已经冻结的认证盒内；锚点、点数、种子、
路由规则和验收条件均未改变。

## 固定验收

1. 全20点严格正/失败与v5最终输出双向一致，FP=FN=0。
2. C1--C4每点必须由预期近似证书实际接受；回退不计为通过。
3. C5/C6必须严格回退；其额外开销单独报告。
4. C1/C2各组含严格NLO前缀的总耗时至少降低20%。
5. C3/C4各组含严格NLO前缀的总耗时必须为正收益；否则即使正确也撤销该局部域。
6. 报告失败阶段漂移和双方positive后的SNR幅值误差，但二者不替代FP/FN门槛。
7. 结果只适用于上述预定义目标分层，不能代表完整BSMPT参数空间。
