# Region-v6 type-I安全收缩报告

日期：2026-09-07

## 原因

冻结的type-I专项E3共18点，正/失败主契约为TP=8、TN=10、FP=FN=0，总降时46.47%。
但是broad-C row3正点域T3的两点均出现严格/近似相变历史schema不一致。两边仍有正SNR，
故不构成FP/FN；但它说明该域不适合继续作为“安全近似”区域。

## v6改动

region-v6只增加一个type-I输入级严格区域：broad-C row3锚点及原`1.1e-4`经验半径。
该规则放在近似首遍之前，避免先付近似成本再回退。其余v5规则和认证半径不变。

为保持旧入口行为，`run_calcgw_approx_region_v3_multitemp.sh`只把exact-anchor路径改成
“环境变量存在时采用调用者值，否则仍用v3默认值”；因此v3/v4/v5默认行为不变。

## 复验

T3两点均打印
`direct exact: known_exact_region:type1_c3_transition_schema_mismatch_strict`。v6与严格输出
除runtime外逐字段完全相同，相变历史schema不一致由v5的2/2降为0/2，FP=FN=0。

当前推荐type-I入口：

```bash
bsmpt_speed_lab/run_calcgw_approx_region_v6_type1_safe.sh \
  --model=r2hdm --input=input.tsv --output=output.tsv \
  --firstline=2 --lastline=2
```

累计严格配对证据258点：TP=96、TN=162、FP=FN=0。新增18点全部为type-I；该累计值
仍包含相关局部扰动，不能外推成完整type-I参数空间保证。
