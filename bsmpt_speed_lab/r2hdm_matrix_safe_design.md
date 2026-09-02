# CalcGW 矩阵优化候选：只组装最终保留的 gauge 上三角

这是一份**不直接应用**的最小补丁建议，目标文件是隔离测试版
`upstream/src/models/ClassPotentialOrigin.cpp`。它只影响 `CalcGW` 使用到的
`GaugeMassesSquared` 路径，不触碰原项目、`install`、数据文件或 thdmtools。

## 为什么选择这一项作为第一步

`GaugeMassesSquared` 当前对 `a,b=0..NGauge-1` 全部组装质量矩阵，随后马上执行：

```cpp
for (a = 1; a < NGauge; ++a)
  for (b = 0; b < a; ++b)
    MassMatrix(a, b) = MassMatrix(b, a);
```

因此每个下三角元素的前一轮求和都会被丢弃。候选补丁在环境变量
`BSMPT_SAFE_GAUGE_UPPER_ONLY=1` 时只计算 `b >= a`，而保留原来的镜像循环。
默认不开启，方便与现有控制版逐点对照。

这不是假设 CP 或电荷守恒，也不是把路径限制在中性场方向：它利用的是当前
基类实现本来就把下三角覆盖成上三角这一事实。对于任意派生模型、CB/CP 路径、
`Temp=0`、有限温度以及 `diff=0, diff>0, diff=-1`，都只减少被丢弃的工作。

## 数值安全性

启用后，每个最终保留的 `MassMatrix(a,b)` 的初始化和 `i,j` 加法顺序完全不变；
下三角仍逐元素复制上三角。因此在 IEEE 双精度下，正常数值结果应当逐 bit
一致，而不是只有物理量近似一致。尤其没有把求和改成 Eigen 矩阵乘法、没有
改变本征值求解器，也没有重排 `V1Loop` 中 boson/fermion 的加法顺序。

唯一需要显式记录的边界是：如果未来代码移除下三角镜像，或者某个派生模型依赖
一个故意非对称的 gauge tensor，那么这个补丁不再适用；在当前基类契约下，后者
不会进入本征值计算，因为旧代码也会覆盖下三角。

## 复杂度和预期收益

矩阵元素组装从 `NGauge^2 * NHiggs^2` 降为
`NGauge*(NGauge+1)/2 * NHiggs^2`。R2HDM 当前 `NGauge=4, NHiggs=8`，这一段
的乘加次数理论上从 1024 降到 640（减少 37.5%）。但它只占一次 `VEff` 的一部分，
而 CalcGW 的主要耗时在 bounce 的大量势能/导数调用，所以总 wall time 很可能
只有小幅下降；必须实测，不应预先宣称大加速。

## 验证方法

在同一个隔离构建中，用完全相同的单点输入分别运行：

```bash
./CalcGW benchmark_input.tsv control.tsv
BSMPT_SAFE_GAUGE_UPPER_ONLY=1 ./CalcGW benchmark_input.tsv upper.tsv
```

实际参数和输出文件名按现有 `bsmpt_speed_lab` 驱动脚本调整。验证顺序：

1. `mH<125` 和 `mH>125` 各至少一个已计算 classified 点；
2. 用已有 `compare_outputs.py --rtol 0 --atol 0` 比较所有非 runtime 字段；
3. 对一个高 SNR 点和一个 SNR 极小点重复运行，检查跃迁数、`Tn/Tp`、`alpha`、
   `beta/H` 和 SNR；
4. 记录 `BSMPT_PROFILE_PHASES=1`，确认只减少矩阵组装时间，且没有改变失败/成功
   状态；
5. 若所有输出逐 bit 一致，再进行至少三次 wall-time 中位数 benchmark。

## 关于“活动非零场索引”与 T=0/有限 T 复用

R2HDM 的 bounce 降维向量通常只写入 `VevOrder={2,4,6,7}`，因此确实存在大量
严格为零的 full-VEV 分量。下一项可以在**另一个独立开关**下跳过
`v[i]==0 || v[j]==0` 的乘积项，并保留原 `i,j` 顺序；不过它会删除对有符号零的
无效加法，不能像本补丁一样直接保证逐 bit 相同，必须先做全量数值回归。因此这里
不把它和第一项混在一起。

同理，`V1Loop` 每次需要有限温度和 T=0 的不同质量谱。把场依赖矩阵缓存后再加
Debye 项理论上更有收益，但当前 API 返回的是已经对角化后的谱，且 `const` 方法、
温度和 `diff` 组合很多；若没有明确的 cache 生命周期和参数失效规则，容易把不同
温度或不同导数混用。应先单独设计 `FieldMassMatrix(v)`/`AddDebye(T)` 接口，确认
所有派生模型都满足分解，再做第二阶段实验。

## 结论

`r2hdm_matrix_safe_proposal.patch` 是当前最安全、最小、可回退的矩阵级候选：默认
零行为变化，启用后严格只省略旧代码最终会覆盖的下三角计算。若 benchmark 收益
过小，应把主要精力转向带独立开关的活动索引或显式 field-matrix cache，而不是
盲目替换 Eigen 本征值算法。
