# BSMPT CalcGW 代码级优化完整报告

最后更新：2026-09-01

## 1. 目标、范围和结论

本工作的目标是缩短 BSMPT `CalcGW` 对单个 R2HDM 点的计算时间，以便把
`neither` 区域中的大量点送入完整 GW 计算。工作只涉及 CalcGW/BSMPT 内部
算法和实现，不涉及 thdmtools，也不把点级并行计入加速结果。

所有源码、构建和输出均位于 `bsmpt_speed_lab`。原程序
`install/bin/CalcGW` 没有被覆盖、替换或重新链接，系统环境和现有软件连接也
没有被修改。

当前推荐结果如下：

| 分支 | 原安装程序 | 当前严格版本 | 时间减少 | 吞吐提升 |
|---|---:|---:|---:|---:|
| `mH < 125` 主参考点 | 203.37 s | 95.490 s | 53.0% | 2.13x |
| `mH > 125` 主参考点 | 88.66 s | 29.900 s | 66.3% | 2.97x |
| 4 个 classified 独立点 | 248.635 s | 126.270 s | 49.2% | 1.97x |
| `mH < 125`、SNR≈150 | 140.172 s | 79.463 s | 43.3% | 1.76x |

`mH<125` 与 `mH>125` 的运行行为确实不同。因此提供两个明确入口，调用者
必须事先按物理分支选择，程序不会自动判断：

- `mH<125`：`run_calcgw_exact_fast.sh`
- `mH>125`：`run_calcgw_exact_fast_mh_gt_125.sh`

除特别标记为“近似”的实验外，严格版本没有减少 threshold 的 1000 点网格、
没有降低有限差分阶数、没有更改积分容差，也没有把完整费米子本征问题降维。

## 2. 源码、构建和数据边界

| 内容 | 路径 |
|---|---|
| 隔离优化源码 | `bsmpt_speed_lab/upstream` |
| 当前正式测试二进制 | `bsmpt_speed_lab/build-conda-control/bin/CalcGW` |
| 原始安装二进制 | `install/bin/CalcGW` |
| 通用严格 wrapper | `bsmpt_speed_lab/run_calcgw_exact_fast.sh` |
| `mH>125` 严格 wrapper | `bsmpt_speed_lab/run_calcgw_exact_fast_mh_gt_125.sh` |
| 输出比较工具 | `bsmpt_speed_lab/compare_outputs.py` |

隔离源码来自 BSMPT 3.3.1 Git 对象
`04cb17d1233522f3c423cbd957a8922be037241e`。`upstream` 本身没有 `.git`
目录，因此不会修改父项目的 Git 历史。

`build-native*`、`build-lto`、`build-pgo-gcc`、`build-nostack-*` 是历史实验
目录，不应替代 `build-conda-control` 作为正式结果来源。

CalcGW 输入必须只有前 8 个模型参数列：

```text
yuktype L1 L2 L3 L4 L5 m12sq tbeta
```

CalcGW 的完整输出不能直接再次作为输入。若需要从旧输出恢复输入，必须先取
前 8 列：

```bash
cut -f1-8 previous_output.tsv > clean_input.tsv
```

早期曾误把完整输出作为输入，形成附加列；相应运行不用于本报告的准确性结论。

## 3. 正确性判据

严格验证忽略且仅忽略 `runtime` 列，对状态字符串和所有其他浮点字段使用零
容差逐字段比较：

```bash
python3 bsmpt_speed_lab/compare_outputs.py \
  --rtol 0 --atol 0 \
  reference.tsv candidate.tsv
```

报告中的“逐位一致”表示上述命令通过，而不是仅比较 SNR 或使用默认容差。

需要区分两种基准：

1. 原安装程序用于衡量实际加速比。
2. 修正后的 dense control 用于严格物理结果比较。

上游路径曾在赋值 `Initial_lmin` 之前使用它，属于未定义行为。修复后结果稳定，
但不要求复现旧二进制偶然产生的最后若干位。修复后的重复运行在两个质量分支
均逐位稳定。此后所有优化都与修正后的 control 比较。

代表性验证文件：

| 用途 | 文件 |
|---|---|
| `mH<125` 输入 | `benchmark_input.tsv` |
| `mH>125` 输入 | `benchmark_input_mh_gt_125.tsv` |
| 四点输入 | `classified_validation_4.tsv` |
| 四点修正基准 | `classified_validation_4_exact_baseline.tsv` |
| 四点严格优化 | `classified_validation_4_combined_derivatives.tsv` |
| SNR≈150 干净输入 | `snr150_input_mh_lt_125.tsv` |
| SNR≈150 修正基准 | `dense_corrected_validation_high_snr_mh_lt_125.tsv` |
| SNR≈150 严格优化 | `full_exact_validation_snr150_clean_mh_lt_125.tsv` |
| 最新 `mH>125` 专用输出 | `mh_gt_specialized_latest.tsv` |

四点零容差复现示例：

```bash
python3 bsmpt_speed_lab/compare_outputs.py \
  --rtol 0 --atol 0 \
  bsmpt_speed_lab/classified_validation_4_exact_baseline.tsv \
  bsmpt_speed_lab/classified_validation_4_combined_derivatives.tsv
```

预期结果：`PASS: 4 rows agree (runtime ignored)`。

## 4. 当前正式严格配置

### 4.1 通用配置：主要用于 `mH<125`

`run_calcgw_exact_fast.sh` 启用：

```text
BSMPT_USE_PATH_GEOMETRY_JET=1
BSMPT_USE_SPLINE_INTERVAL_HINT=1
BSMPT_USE_COMBINED_NUMERICAL_DERIVATIVES=1
BSMPT_USE_R2HDM_HIGGS_PAIR=1
BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF=1
BSMPT_USE_R2HDM_HIGGS_INDEX_CACHE=1
BSMPT_USE_R2HDM_GAUGE_PAIR=1
BSMPT_USE_R2HDM_GAUGE_INDEX_CACHE=1
BSMPT_USE_R2HDM_QUARK_INDEX_CACHE=1
BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE=1
BSMPT_USE_R2HDM_QUARK_FIXED12_DIFF0=1
BSMPT_USE_R2HDM_LEPTON_FIXED9_DIFF0=1
BSMPT_USE_V1LOOP_THERMAL_CONTEXT=1
BSMPT_USE_R2HDM_VTREE_INDEX_CACHE=1
BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE=1
BSMPT_USE_R2HDM_COUNTERTERM_FLAT_CACHE=1
BSMPT_USE_FERMION_LOW4_EXACT=1
BSMPT_USE_V1LOOP_MASS_VECTOR_RESERVE=1
BSMPT_USE_CACHED_PROFILE_GATE=1
```

调用方法：

```bash
./bsmpt_speed_lab/run_calcgw_exact_fast.sh \
  --model=r2hdm \
  --input=mh_lt_125_points.tsv \
  --output=gw_results.tsv \
  --firstline=2 --lastline=100
```

### 4.2 `mH>125` 兼容入口

`run_calcgw_exact_fast_mh_gt_125.sh` 现直接调用通用配置，不再增加分支专用开关。

调用方法：

```bash
./bsmpt_speed_lab/run_calcgw_exact_fast_mh_gt_125.sh \
  --model=r2hdm \
  --input=mh_gt_125_points.tsv \
  --output=gw_results.tsv \
  --firstline=2 --lastline=100
```

这两个 fixed-size 路径仍分别求完整的复 12×12 quark 谱和复 9×9 lepton
谱；没有使用 6×6 分块近似，也没有把数学上简并的本征值复制两次。因此 Eigen
求解、排序和近零截断结果与动态矩阵路径逐位一致。

删除动态中间矩阵并直接构造 fixed-size fermion 矩阵后，这两个完整求解路径在
`mH<125` 也出现稳定收益，因此已移入通用 wrapper。当前验证中主点 fixed9
路径为 95.531 s，加入 fixed12 后为 93.511 s；四点集合从 126.270 s 降至
121.957 s。两个主点、四个 classified 点和 SNR≈150 点的所有非 runtime
字段均逐位一致。flat CounterTerm 也已进入通用 wrapper。

## 5. 已接受优化的详细说明

### 5.1 删除必然被覆盖的 bounce 栅格

`BounceActionInt::SetPath` 原先构造 1001 点 `dV/dl` 栅格，但
`Solve1DBounce` 在第一次使用前总会重新生成它。删除这段死工作后输出逐位一致，
早期主点约从 201.9 s 降至 188.7 s。

### 5.2 修复 `Initial_lmin` 和复用真空 Hessian

先设置 `Initial_lmin`，再进行 exact threshold 扫描，消除未定义行为。收敛后的
真空 Hessian 在同一路径上是常数，因此只计算一次。这个修复既改善稳定性，也是
早期最大的一项代码级削减：当时 `mH<125` 约降至 145.47–146.17 s，
`mH>125` 约降至 55.93–60.71 s。

### 5.3 合并数值梯度和 Hessian

原算法分别计算四点梯度与完整 Hessian。合并路径保留全部原 stencil 和相同
`eps`，仅复用二者共有的 `phi ± 2 eps e_i` 势函数点。没有使用二点 central
gradient，也没有删除混合 Hessian 元素。这一项最终把势函数调用数显著压低，是
当前 95.490 s/30.535 s 阶段的关键优化。

### 5.4 路径几何 jet 与 spline interval hint

路径几何 jet 在一次 spline 反演中同时得到 `phi`、`dphi/dl`、
`d2phi/dl2`。代表性配对：

- `mH>125`：61.721 → 57.601 s；
- `mH<125`：157.033 → 142.068 s。

spline interval hint 为连续积分查询保留上一次区间；查询越界或方向改变时仍会
回退查找，不跨路径保存状态。

### 5.5 Simpson 端点复用

最终势能 action 的相邻 Simpson 区间共享端点。新实现复用已经计算的端点值，
不改变采样坐标或加法次序；收益较小但严格安全。

### 5.6 Higgs 和 Gauge 的 T/0 成对矩阵构造

在 Arnold–Espinosa 路径中，`M(v,T)` 与 `M(v,0)` 的场依赖部分相同，
仅 Debye 热项不同。优化只共享矩阵构造，仍执行两个独立 Eigen 求解。

- Higgs `diff=0` 代表性收益约 5.4%；
- Higgs `diff=1..8` 在两个质量分支和四个 classified 点均逐位一致；
- derivative pair 的增量收益约为 `mH>125` 1.1%、`mH<125` 6.5%、
  四点总计 2.3%；
- Gauge pair 遵守相同原则，不共享谱结果。

### 5.7 有序非零索引缓存

Higgs、Gauge、Quark、Lepton、VTree 和 CounterTerm 的高阶张量包含大量恒零
项。缓存只记录每个模型实例中的有序非零索引，不缓存任何随 VEV、温度或参数点
改变的数值。非零项仍按原嵌套循环顺序累加，因而避免了重新排序带来的舍入差异。

缓存会在模型重新初始化或参数改变后重建，不跨参数点共享。

### 5.8 V1Loop thermal context 与热系数连续存储

一次 `V1Loop` 内复用 `T²`、`T³`、`T⁴` 和 `2π²`。热展开系数从
`std::map<int,double>` 改为连续 vector，并把已知范围内的访问内联为直接索引。
所有 `pow/log/sqrt` 表达式和求和顺序保持不变。细分计时约减少 0.18 s，
总墙钟收益较小。

### 5.9 跳过禁用的 `VTreeSimplified`

R2HDM 明确设置 `UseVTreeSimplified=false`。旧代码仍虚调用一个恒返回 0、随后
丢弃结果的函数；新路径直接进入显式树势。输出逐位一致，单项收益小于系统波动。

### 5.10 flat CounterTerm

把已有 CounterTerm 层级索引展开成保持原 `i→j→k→l` 顺序的 flat stream，
不预乘系数、不合并项。两轮 `mH>125` 测试平均约快 2.5%，逐位一致。直接
fixed fermion 路径加入后重新测试，`mH<125` 主点从 93.511 s 降至 92.476 s，
高 SNR 点从 69.718 s 降至 68.679 s，四点集合从 121.957 s 降至 120.622 s；
因此现已在通用 wrapper 中启用。

### 5.11 `mH>125` 完整固定尺寸 fermion Eigen

动态矩阵改为固定编译期尺寸后仍执行完整矩阵乘法和完整
`SelfAdjointEigenSolver`：

- quark fixed 12×12 单独：31.584 s；
- lepton fixed 9×9 单独：31.616 s；
- 同版本动态对照：32.336 s；
- 两者组合：30.679 s；
- 清除 profiling 判断后，同版本配对 31.48 → 30.50 s；
- 与 flat CounterTerm 组合：29.900 s。

所有非 runtime 字段的 SHA-256/零容差逐字段比较一致。

## 6. 性能剖析和剩余瓶颈

`mH>125` 严格配置的阶段剖析：

| 阶段 | 时间 |
|---|---:|
| vacuum/phase tracing | 8.316 s |
| bounce construction | 26.894 s |
| nucleation/percolation/completion/GW 后处理 | <0.14 s |

bounce 内部：

| 子阶段 | 调用数 | 时间 |
|---|---:|---:|
| `Solve1DBounce` | 24 | 21.523 s |
| path check | 12 | 1.352 s |
| action integral | 12 | 0.664 s |
| path deformation | 12 | 0.442 s |

VEff 组件剖析：

| 组件 | 聚合时间 |
|---|---:|
| quark spectra | 10.93 s |
| lepton spectra | 5.83 s |
| fermion thermal functions | 2.69 s |
| CounterTerm | 2.61 s |
| boson thermal functions | 1.34 s |
| VTree | 0.78 s |

因此剩余严格瓶颈主要是 dense threshold scan 反复触发的完整 fermion 谱和势
函数计算。最终 nucleation/GW 后处理本身不是值得优先优化的部分。

profiling 开关只用于诊断，不能用于正式 benchmark：

```text
BSMPT_CALCGW_PROFILE
BSMPT_CALCGW_PROFILE_VEFF
BSMPT_PROFILE_PHASES
BSMPT_PROFILE_ACTION
```

细粒度 RAII 计时器曾在关闭状态下仍产生数千万次判断，现已从正式热路径移除。

## 7. 严格但无收益或被更好方案替代的实验

| 实验 | 代表结果 | 结论 |
|---|---:|---|
| Quark eigensolver workspace | 32.72 vs 32.43 s | 慢约0.9%，代码已移除 |
| RK5 vector workspace | 32.684 vs 31.623 s | 变慢，代码已移除 |
| RK5 fixed array | 约44.99 s | 明显变慢 |
| Higgs/Gauge solver 对象复用 | 31.693 vs 31.608 s | 无收益，代码已移除 |
| fixed Higgs/Gauge eigensolver | 41.995 s | 变慢 |
| pair stack matrices | 37.68 s | 变慢 |
| pair output vector reserve | 32.657 s | 无收益，代码已移除 |
| V1Loop hash contribution cache | 33.086 s | 查找开销超过复用收益 |
| V1Loop linear contribution cache | 35.114 s | 变慢 |
| Higgs derivative upper-only | 52.121 s | 变慢 |
| 跳过 Gauge 零项 | 50.337 s | 变慢 |
| logistic 饱和分支提前退出 | 32.870 s | 无收益，代码已移除 |
| const-reference callback | 31.601 s | 无可靠收益 |
| active Hessian dimensions | 4 个方向全部非零 | 无元素可跳过 |
| 关闭热点 stack protector | 38.645 s | 逐位一致但明显变慢 |
| LTO | 收益不稳定 | 曾产生数值变化，拒绝 |
| native vectorization | 不适用 | 发生 heap corruption，禁止 |

`BSMPT_SKIP_ZERO_R2HDM_HIGGS_TERMS` 本身可以严格工作，但已经被有序 Higgs
index cache 覆盖，不再作为推荐入口。

## 8. 近似模式：可筛选，但不能作为最终 GW 结果

### 8.1 二点 central gradient

| 参考点 | 原时间 | 快速时间 | 主要影响 |
|---|---:|---:|---|
| `mH<125`，SNR≈11.55 | 181.56 s | 129.35 s | SNR相对变化约−1.07e−5 |
| `mH<125`，SNR≈150 | 304.54 s | 206.78 s | SNR相对变化约−2.93e−5 |
| `mH>125`，历史SNR≈54.69点 | 619.50 s | 323.66 s | SNR相对变化约−1.08e−3 |

主跃迁大体一致，但第三个点的一个极弱次级跃迁改变状态。因此 central2 不能用于
严格探查所有 `SNR>0` 边界。另需注意：该 `mH>125` 点在后来修正后的严格导数
下实际总 SNR 约为 `2.304e-7`，历史约54.69的值来自近似导数路径，不能作为
当前严格基准。

### 8.2 Adaptive threshold

threshold 样本数可减少约 68%，但独立 SNR≈150 点出现：

- beta/H 约 3% 偏差；
- SNR 约 1.7% 偏差。

因此不能加入严格模式。历史默认 1000 点本身也存在离散误差；改变它属于收敛性
研究，而不是严格等价的代码优化。

### 8.3 解析 lepton 谱

R2HDM lepton 9×9 矩阵数学上分成三个 3×3 块，可以解析得到每代
`0, |a|², |a|²+|b|²`。但解析公式与完整 Eigen 的舍入不同：

| 点 | 数值版 | 解析版 | 误差 |
|---|---:|---:|---|
| `mH<125`，SNR≈150 | 79.46 s | 68.50 s | SNR −2.74%，beta/H −0.757% |
| `mH>125` 测试点 | 141.42 s | 115.77 s | SNR +0.086%，beta/H +0.0165% |

状态字段一致，但 `mH<125` 偏差不可接受。开关
`BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0` 仅保留研究用途。

### 8.4 Quark 12×12 → 6×6

R2HDM quark 矩阵具有严格手性块结构，数学谱由一个 6×6 奇异值问题重复两次。
然而完整 12×12 Eigen 会分别数值计算两个理论简并块；降维后复制本征值改变末位
舍入和后续 V1Loop 累加。实测不仅更慢，还使 beta/H 偏移约 0.31%，因此拒绝。

### 8.5 内部多线程

约从 203.3 s 降至 196.6 s，但总 SNR 相对变化约 `2.65e-5`、beta/H 约
`3.06e-4`。此外用户的扫描本来已经并行，因此它既不是本项目需要的代码级核心
优化，也不适合边界敏感分类。

## 9. PGO 和编译实验

独立 GCC PGO 曾得到：

| 分支 | control | PGO | 结果 |
|---|---:|---:|---|
| `mH>125` | 32.266、32.100 s | 30.407、31.017 s | 平均约4.6%，逐位一致 |
| `mH<125` | 99.146 s | 96.274 s | 约2.9%，逐位一致 |

随后用与 control 完全匹配的 vectorization flags 重建，收益降至约 2%，说明 PGO
对训练点、构建 flags 和机器状态敏感。目前保留 `build-pgo-gcc` 作为独立实验，
不把它计入 29.900 s 的正式代码级结果。

## 10. 当前保留但默认关闭的主要研究开关

以下开关存在于源码，但不在正式 wrapper 中：

- 导数/积分：`BSMPT_USE_CENTRAL2_GRADIENT`、
  `BSMPT_USE_CENTRAL2_DVDL`、`BSMPT_USE_DIRECTIONAL_DVDL`、
  `BSMPT_USE_ANALYTIC_GRADIENT`、`BSMPT_ADAPTIVE_THRESHOLD`、
  `BSMPT_DENSE_THRESHOLD_STEPS`、`BSMPT_BOUNCE_RASTER_INTERVALS`；
- bounce 实验：`BSMPT_USE_RK5_FIXED_STORAGE`、
  `BSMPT_USE_BESSEL_RECURRENCE`、`BSMPT_REUSE_RASTER_IN_LIMIT_SCAN`、
  `BSMPT_USE_ACTIVE_COMBINED_HESSIAN`；
- 矩阵实验：`BSMPT_USE_ANALYTIC_GAUGE_MASSES`、
  `BSMPT_USE_R2HDM_FIXED_EIGENSOLVER`、
  `BSMPT_USE_R2HDM_PAIR_STACK_MATRICES`、
  `BSMPT_USE_R2HDM_QUARK_BLOCK6`、
  `BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0`；
- cache 实验：`BSMPT_USE_V1LOOP_MASS_CONTRIBUTION_CACHE`、
  `BSMPT_USE_V1LOOP_LINEAR_MASS_CONTRIBUTION_CACHE`。

设置这些变量不代表结果经过严格认证。除正式 wrapper 中明确列出的开关外，都应
重新进行两个质量分支、classified 点和高 SNR 点验证。

## 11. 推荐工作流程

1. 在 classified 数据中先按 `mH<125` 和 `mH>125` 分开生成只含 8 个输入列
   的文件。
2. 两个质量分支均可使用通用严格 wrapper。
3. `mH>125` 专用脚本现为兼容入口，调用同一通用配置。
4. 每次源码或编译 flags 改变后，至少验证：两个主点、四个独立 classified 点、
   一个高 SNR `mH<125` 点。
5. 使用 `compare_outputs.py --rtol 0 --atol 0` 检查严格模式；近似筛选模式则
   另外报告状态、温度、alpha、beta/H 和各探测器 SNR 的相对误差。
6. 不把 profiling 运行时间与正式运行混用。

## 12. 安全注意事项

- 不要在未隔离 `CONAN_HOME` 时运行上游 `Build.py`；它可能删除并重建现有
  `~/.conan2/profiles/BSMPT`。
- 不使用发生过 heap corruption 的 native-vectorization 构建。
- 不写入或替换 `install/bin/CalcGW`。
- 不把实验 wrapper 或环境变量加入正式扫描脚本，除非完成跨分支严格验证。
- 所有输出应写入 `bsmpt_speed_lab` 或用户明确指定的数据目录。

## 13. 后续最有价值的研究方向

### 13.1 2026-09-01 后续严格优化：fixed fermion 直接构造与 flat CT

细粒度 profile 显示，在 `mH>125` 主参考点约 227 万次 VEff 调用中，quark
和 lepton 质量路径分别累计约 9.77 s 和 5.19 s。进一步检查发现：即使正式
wrapper 已选择固定 12x12 quark 和 9x9 lepton 路径，函数入口仍会无条件分配
随后不用的动态 Eigen 矩阵。现改为仅在动态 fallback 被实际选择时才分配，并
直接按原索引顺序构造固定矩阵；矩阵元素、完整 Eigen 求解器、截断条件和浮点
累加顺序均未改变。

该修改在两个主参考点、四个 classified 独立点和 `mH<125` SNR≈150 点上均
通过 `--rtol 0 --atol 0` 比较，共 7 个点的所有非 runtime 字段逐位一致。并发
fixed9 在 `mH<125` 主点为 95.531 s，加入 fixed12 后为 93.511 s；高 SNR
点对应为 72.491 s 和 69.718 s。四点集合在 fixed9+fixed12 下为 121.957 s，
再加入 flat CounterTerm 后为 120.622 s，对比此前 126.270 s。flat CT 在
`mH<125` 主点和高 SNR 点也分别进一步降至 92.476 s 和 68.679 s，因此三个
开关均已移入通用 wrapper。并发负载会影响绝对时间，候选间的同阶段比较用于
判断方向，不覆盖报告前部较稳定的历史基准。

本轮还复测并拒绝了四项逐位一致但无收益的候选：Higgs 单动态工作矩阵、Higgs/
gauge 固定尺寸求解器、显式缓存 thermal-mass `Temp^2`、精确零质量 fermion
thermal 单值缓存。失败候选代码已从正式路径移除或保持既有实验开关关闭。

当前没有尚未筛选且证据充分的低风险严格候选。剩余研究空间集中在：

1. Eigen 上游或编译器层面对完整 12×12/9×9/8×8 数值路径的改进；
2. 能证明不改变 `pow/log/sqrt` 和累加顺序的 thermal function 专门化；
3. dense threshold scan 中跨 stencil 的模型内部公共量复用；
4. 分别针对 `mH<125` 和 `mH>125` 建立更大的 paired benchmark 集，而不是
   假设一个开关在两侧都有相同收益。

这些方向目前都需要新的算法证据或更大的 paired 数据集，不能作为现成的严格
优化继续启用。任何通过减少 threshold 点数、降低有限差分阶数或复制解析简并本征值获得的速度，
都应标记为近似筛选，不计入严格 CalcGW 加速。

### 13.2 固定四项 fermion 低温展开

`JfermionInterpolatedLow` 的正式调用阶数固定为 4，但编译器仍保留运行时循环。
实验室版本按原来的 `l=2,3,4` 顺序显式写出三项，保留相同的 `pow` 调用、系数
和加法顺序。两个主点、四个 classified 点以及 SNR≈150 点均通过非 runtime
字段 `--rtol 0 --atol 0` 比较。

同负载 paired 测试中，`mH>125` 两轮分别为 30.594 vs 31.214 s、31.581 vs
31.679 s；高 SNR 点为 73.395 vs 74.449 s。收益约为 0.3%--2.0%，明显小于
此前结构性优化，但三组方向一致，因此 `BSMPT_USE_FERMION_LOW4_EXACT=1` 已加入
实验室正式 wrapper。不同时间或并发状态下的绝对 runtime 不用于计算该增益。

### 13.3 V1Loop 成对谱结果容量预留

正式 R2HDM Higgs pair 路径的两个 8 元素结果，以及 gauge pair 路径的两个
4 元素结果，原先在逐个 `push_back` 时发生多次扩容。候选仅按 `NHiggs` 或
`NGauge` 预留容量；矩阵构造、完整 Eigen 求解、阈值处理与元素写入顺序均不变。

两个主点 paired 测试分别为 30.527 vs 30.821 s（`mH>125`）和 95.433 vs
96.280 s（`mH<125`），收益约 0.95% 和 0.88%。两个主点、四个 classified
点和 SNR≈150 点共 7 点均通过非 runtime 字段零容差验证，因此
`BSMPT_USE_V1LOOP_MASS_VECTOR_RESERVE=1` 已加入实验室正式 wrapper。

### 13.4 关闭 profiling 时缓存计时 gate

VEff 细粒度计时开关本身是进程级静态值，但原先每个 boson/fermion contribution
仍经非内联函数重复查询。候选在 V1Loop 内缓存同一个 gate，并把它传给计时器；
开启 profiling 时仍读取时钟并记录所有原指标，关闭时只省去重复查询。

`mH>125` paired 为 30.599 vs 31.144 s，`mH<125` 为 96.703 vs
97.004 s，分别约快 1.75% 和 0.31%。共 7 点非 runtime 字段零容差一致，故
`BSMPT_USE_CACHED_PROFILE_GATE=1` 已加入实验室正式 wrapper。该收益较小，仍可能
受调度噪声影响，不与早期结构性加速混算。

### 13.5 当前 Eigen/V1Loop 深层审计结论

- quark 与 lepton 正式路径已经直接构造固定 12×12/9×9 矩阵并运行完整
  `SelfAdjointEigenSolver`；动态 fallback 矩阵不会在该路径分配。
- Higgs/gauge 的固定尺寸求解器、单个动态 solver workspace、pair stack matrix
  均已做 paired 测试但无收益，继续复用 Eigen 对象没有可信证据。
- V1Loop 的 quark/lepton `thread_local` 输出缓冲复用虽逐位一致但略慢，未启用；
  当前新增的容量预留只作用于正式 Higgs/gauge pair 的短结果 vector。
- contribution cache 的四个容器在缓存关闭时不发生堆分配；此前 hash/linear
  cache 候选也无稳定收益，因此不再为它引入指针或 optional 分支。
- 任何 6×6 分块、解析简并复制、解析 lepton 谱或降低有限差分/threshold 点数
  的方案都会改变完整 Eigen/数值路径，继续保留为近似研究项，不计入严格配置。

因此，在保持完整谱求解、原浮点运算顺序和零容差门槛的约束下，本轮能由现有
证据支持的 V1Loop 分配与 Eigen 调用优化已经筛完。进一步显著收益需要新的
算法等价证明、Eigen 上游实现变化或更大 paired 数据集，而不是重复已有候选。

### 13.6 三个新增独立 GW 点与 thermal 消融

从此前未用于 CalcGW 接受门槛的 R2HDM 单元测试参数中筛出三个完整 GW 点，
单独分类运行约为 43.7--52.0 s，LISA SNR 覆盖约 `1.6e-26` 到 `2.8e-7`。
输入保存在 `extra_validation_active_3.tsv`。

当前 wrapper 对比关闭最近三项微优化的 exact control，逐点分别快 2.26%、
3.00%、3.23%，合计 149.445 vs 153.847 s（2.86%）；所有非 runtime 字段
零容差一致。只关闭 fermion low4 展开时，三点合计 148.408 vs 151.627 s，
当前配置快 2.12%；同时关闭结果容量预留与 profiling gate 时，三点合计
136.773 vs 138.517 s，当前配置快 1.26%。两次消融中三个点均同方向。

另测试了固定三阶 boson 低温展开：三个新增点合计快 0.83%，两个质量主点各约
快 0.46%--0.48%，classified 四点也逐位一致；但 SNR≈150 点出现 53 个非
runtime 字段差异（例如 `beta/H` 从 59.6391690163 变为 59.6335080005）。
因此该候选未通过严格门槛，源码和 wrapper 均已撤销。

### 13.7 第二批八个独立点的扩展正确性验证

又从未参与开发或此前验收的 R2HDM 单元测试参数中选择八点，保存在
`extra_validation_candidates_8_v2.tsv`。其中七点完整得到 GW，单点运行约
38--72 s，LISA SNR 覆盖约 `2.1e-28` 到 `1.7e-8`；另一个点通过 NLO 和
tracing，但 GW 状态为 `nan`，作为边界状态样本保留。

当前 wrapper 与 pre-micro exact control 在相同双进程负载下对照，八点所有
非 runtime 字段均通过 `--rtol 0 --atol 0`。七个有效 GW 点合计为
468.903 vs 476.674 s，当前配置快 1.63%；逐点有六个加速、两个受调度噪声或
路径占比影响变慢，因此不声称微优化对每个单点都必然降低 runtime。

随后从中取四个代表点，对比关闭全部 exact-fast opt-in 的 corrected dense
control。四点仍全部逐位一致，耗时分别降低 67.11%、52.74%、51.91%、
51.42%；合计 170.037 vs 392.366 s，降低 56.66%，速度为 2.308×。这组
独立样本确认累计结构性加速并非旧 benchmark 特例。dense control 入口为
`run_calcgw_dense_control.sh`，不会读取或写入原版 `install/`。

### 13.8 NLO 有效性边界压力测试

以一个完整 GW 有效点和一个 `no_nlo_stability` 点为端点，在七维输入参数空间
做线性插值。粗扫发现状态并非单调：`t=0.68` 为完整 success，`t=0.69--0.70`
为 NLO 无效，`t=0.71` 又为完整 success，`t=0.72` 无效，`t=0.74` 有效，
`t=0.76` 再次无效。这种交替状态比普通随机点更适合检验边界稳定性。

选择 `t=0.64、0.66、0.68、0.71、0.74` 五个边界内侧有效点，对比当前 wrapper
与关闭全部 exact-fast opt-in 的 dense control。五点所有非 runtime 字段均
逐位一致，并在两侧都保持 NLO=`success`、GW=`success`。单点 runtime 降低
52.03%--52.38%；合计 183.053 vs 383.043 s，降低 52.21%，即 2.093×。

另对紧邻的 `t=0.69、0.70、0.72、0.76` 四个边界外侧点做相同对照；两侧均
严格得到 `no_nlo_stability / tracing=not_set`，所有非 runtime 字段逐位一致。
因此优化既没有改变边界内侧结果，也没有错误放行边界外侧点。最终输入分别保存
为 `nlo_boundary_valid_5.tsv` 与 `nlo_boundary_invalid_4.tsv`。

### 13.9 广域分层 A/B 组

为避免验证集中在同一条插值线上，又从未使用的 R2HDM 参数中建立两个各十点的
广域组。A 组混合低/中 `m12²`、较高 `tanβ`、正负 `L5` 和强耦合；B 组集中
覆盖大 `m12²`、低 `tanβ` 与更极端耦合。两组都直接对比当前 exact-fast 与
关闭全部 opt-in 的 corrected dense control。

A 组 10/10 非 runtime 字段零容差一致，覆盖完整 GW、GW=`nan`、
`no_nlo_stability`、`non_bfb`、以及 tracing success 但 GW=`not_set`。合计
runtime 为 387.929 vs 785.421 s，降低 50.61%，即 2.025×。

B 组同样 10/10 零容差一致，包含两个 `no_nlo_stability`、多个 `non_bfb`、
完整 GW、GW=`nan` 与 GW=`not_set`，并包含 optimized 单点 176 s 和 405 s
的长尾路径。合计为 825.687 vs 1629.099 s，降低 49.32%，即 1.973×。

输入与输出分别以 `stratified_broad_group_a_10*`、
`stratified_broad_group_b_10*` 命名。至此新增广域/边界测试已覆盖快速拒绝、完整
GW 和数百秒长尾，不再仅依赖原来的七点接受门槛。

### 13.10 广域 C 组与 Yukawa 类型覆盖

C 组再加入十个独立 R2HDM 点，重点覆盖 `L3<0`、较大的正负 quartic 以及更宽的
质量分支。10/10 行与 dense control 的所有非 runtime 字段零容差一致；其中五点
完整得到 GW、四点为 GW=`nan`、一点为 `no_nlo_stability`，并包含 optimized
单点约 313 s 的长尾。合计 runtime 为 556.091 vs 1055.075 s，降低 47.29%，
即 1.897x。输入与结果以 `stratified_broad_group_c_10*` 命名。

为避免完整 fermion Eigen 路径只在 Yukawa type 1 下得到验证，又选择三个标量点，
分别改为 Yukawa type 2、3、4，共九行。9/9 行严格一致；六行为
EWSR/GW success，三行为 `non_bfb` 但仍完成 GW，SNR 从约 `1.5e-26` 到
`2.83e-7`。合计 405.369 vs 810.131 s，降低 49.96%，即 1.999x。该组文件为
`stratified_yukawa_types_9*`，直接覆盖了 type 2/3/4 下完整 12x12 quark 与
9x9 lepton Eigen 路径。

随后把已知 SNR 约 150 的标量点改为 Yukawa type 2/3/4。三行全部严格一致，
type 3 保持 SNR=150.868，而 type 2/4 落入约 `5.4e-13` 的弱信号分支，因此同一
标量点同时覆盖强弱信号路径。合计 runtime 为 220.162 vs 430.514 s，降低
48.86%，即 1.955x。文件以 `stratified_high_snr_yukawa_3*` 命名。

### 13.11 CalcGW multistepmode 模式矩阵

在上述强信号 type-3 点上分别运行 `multistepmode=0、1、2、auto`，逐一对比当前
exact-fast wrapper 与 dense control。四种模式均通过零容差比较，非 runtime
字段 4/4 完全一致：

| mode | optimized / s | dense / s | runtime 降低 |
|---|---:|---:|---:|
| 0（single-step） | 70.471 | 142.266 | 50.47% |
| 1（tracing coverage） | 72.175 | 143.941 | 49.86% |
| 2（global-minimum tracing） | 143.627 | 272.854 | 47.36% |
| auto | 71.760 | 143.209 | 49.89% |

四种模式合计 358.033 vs 702.271 s，降低 49.02%。模式 2 的运行时间约为其他
模式两倍，说明测试实际进入了更重的 global-minimum tracing 路径。对应结果以
`multistepmode*_high_snr_*` 命名。本轮累计至少 73 个新行级 A/B 对照通过严格
门槛，状态覆盖 success、`nan`、`not_set`、`non_bfb`、`no_nlo_stability`，
并包含有效性边界、不同 Yukawa 类型、强/弱 SNR 与四种 multistep 模式。
