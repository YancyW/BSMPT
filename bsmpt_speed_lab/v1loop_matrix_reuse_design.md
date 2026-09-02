# V1Loop 中 Higgs/Gauge 的 T 与 T=0 质量矩阵复用设计

本文是对 `upstream` 隔离测试版的只读设计记录；本轮没有修改
`bsmpt_speed_lab/upstream` 中的源码，也没有修改原项目。

## 结论

在一次 `Class_Potential_Origin::V1Loop(v,T,diff)` 调用内，当前代码会把同一
个场点的场依赖质量矩阵分别组装多次：

```text
M_H(T), M_G(T), M_G(0), [M_H(0)]
```

非 Parwani 分支的方括号项会被调用；Parwani 分支通常不需要 `M_H(0)`，但
仍然需要 `M_G(0)`。质量矩阵的温度依赖只有最后的 Debye 项：

```text
M_H(T) = M_H(0) + T^2 DebyeHiggs
M_G(T) = M_G(0) + T^2 DebyeGauge
```

因此可以在一次 `V1Loop` 中只组装一次场依赖部分，复制得到 `T=0` 与 `T`
矩阵，然后分别调用现有 `SelfAdjointEigenSolver`。这不改变本征值算法，也
不降低矩阵维数或积分精度。

预期只减少矩阵组装，不会消除本征值求解：

| V1Loop 分支 | 当前矩阵组装 | 复用后 | 可去掉的场依赖组装 |
|---|---:|---:|---:|
| `diff==0`, 非 Parwani | H(T), G(T), G(0), H(0) | H pair + G pair | 2 次 |
| `diff==0`, Parwani | H(T), G(T), G(0) | H(T) + G pair | 1 次 |
| `diff>0`, 非 Parwani | H(T), H(0), G(T), G(0) | H pair + G pair | 2 次 |
| `diff>0`, Parwani | H(T), G(T), G(0) | H(T) + G pair | 1 次 |
| `diff==-1`, 非 Parwani | H(T), H(0), G(T), G(0) | H pair + G pair | 2 次 |
| `diff==-1`, Parwani | H(T), G(T), G(0) | H(T) + G pair | 1 次 |

这里“去掉”仅指重复的多项式矩阵组装；两套矩阵仍各自送入原有
eigensolver，以避免改变求解器数量和行为。

## 当前代码路径与可利用的不变量

`V1Loop` 先调用：

```cpp
HiggsMassesSquared(v, Temp, diff);
GaugeMassesSquared(v, Temp, diff);
GaugeMassesSquared(v, 0, diff);
QuarkMassesSquared(v, diff);
LeptonMassesSquared(v, diff);
```

非 Parwani 的 `diff==0` 分支另外调用 `HiggsMassesSquared(v, 0)`；非
Parwani 的 `diff>0` 和 `diff==-1` 也需要 Higgs 的 `T=0` 质量。

`HiggsMassMatrix(v,Temp)` 的场依赖部分是

```cpp
L2(i,j)
+ sum_k L3(i,j,k) * v[k]
+ sum_{k,l} 0.5 * L4(i,j,k,l) * v[k] * v[l]
```

温度非零时只在该和的末尾增加 `DebyeHiggs(i,j) * pow(Temp,2)`。
`GaugeMassesSquared` 具有完全相同的结构，场依赖部分为

```cpp
sum_{i,j} 0.5 * G2H2(a,b,i,j) * v[i] * v[j]
```

温度项是 `DebyeGauge(a,b) * pow(Temp,2)`。

在一轮 CalcGW 单点计算中，`Curvature_*` 和 Debye 矩阵在 `set_All`/
`initModel` 后不应改变；`v` 和 `Temp` 是该次 `V1Loop` 的固定输入。因此
这个代数拆分成立。对 `diff>0`，场质量矩阵的导数矩阵也与温度无关；对
`diff==-1`，导数矩阵是 `2*T*Debye`，不能把 T 与 0 的导数矩阵混为一谈。

## 最小 API 改动建议

推荐把缓存限制在 `V1Loop` 的局部对象，不添加跨调用的 mutable 成员缓存。
可在 `Class_Potential_Origin` 的 `protected` 区域增加如下内部类型和辅助
函数（以下为伪码，未实际应用）：

```cpp
struct ThermalMassMatrixPair {
  Eigen::MatrixXd zero;    // M(v,0)，只在 needZero 时有效
  Eigen::MatrixXd thermal; // M(v,T)
  bool hasZero = false;
};

ThermalMassMatrixPair BuildHiggsMassMatrixPair(
    const std::vector<double>& v, double Temp, bool needZero) const;

ThermalMassMatrixPair BuildGaugeMassMatrixPair(
    const std::vector<double>& v, double Temp, bool needZero = true) const;

std::vector<double> EigenvaluesSquared(
    const Eigen::MatrixXd& M, double zeroMass) const;

std::vector<double> EigenvalueDerivatives(
    const Eigen::MatrixXd& M, const Eigen::MatrixXd& dM) const;
```

如果希望改动更小，`EigenvaluesSquared` 和 `EigenvalueDerivatives` 可以不
抽出，直接把目前 `HiggsMassesSquared`/`GaugeMassesSquared` 中的 eigensolver
循环复制到 `V1Loop` 的局部 lambda；但抽出两个小 helper 可避免重复实现，且
不会改变 solver 类型。

另一种更保守的实现是保留现有 `HiggsMassMatrix(v,0)`，只增加一个内部
`AddHiggsDebye(MatrixXd&,Temp)`；Gauge 则将当前矩阵组装循环提取为内部
`BuildGaugeMassMatrix(v,Temp)`。这比把矩阵缓存做成 public API 更安全。

## Pair 组装伪码

### Higgs

```cpp
pair BuildHiggsMassMatrixPair(v, T, needZero) {
  // 与现有 HiggsMassMatrix(v,0) 完全相同的 L2/L3/L4 循环和顺序
  MatrixXd M0 = BuildHiggsFieldMatrix(v);
  MatrixXd MT = M0;

  if (T != 0) {
    for (i = 0; i < NHiggs; ++i)
      for (j = i; j < NHiggs; ++j)
        MT(i,j) += DebyeHiggs[i][j] * std::pow(T, 2);
  }

  // 与当前 HiggsMassMatrix 一样把下三角复制自上三角
  SymmetrizeLower(M0);
  SymmetrizeLower(MT);

  return {needZero ? std::move(M0) : MatrixXd{}, std::move(MT), needZero};
}
```

如果 `needZero==false`，仍然可以在内部先得到 `M0` 再丢弃；更好的实现是
保持 `M0` 作为场依赖工作矩阵，但只在需要时保留返回值。不能直接把
`MT` 当成 `M0` 的引用，因为 `diff` 分支及后续调用需要两份独立矩阵。

### Gauge

```cpp
pair BuildGaugeMassMatrixPair(v, T, needZero) {
  MatrixXd M0(NGauge, NGauge);

  if (!SafeGaugeUpperOnly) {
    // 必须按现有代码对每个 (a,b) 独立求和；不能擅自只算上三角。
    for (a = 0; a < NGauge; ++a)
      for (b = 0; b < NGauge; ++b)
        M0(a,b) = GaugeFieldSum(a,b,v);
  } else {
    for (a = 0; a < NGauge; ++a)
      for (b = a; b < NGauge; ++b)
        M0(a,b) = GaugeFieldSum(a,b,v);
    SymmetrizeLower(M0);
  }

  MatrixXd MT = M0;
  if (T != 0) {
    if (!SafeGaugeUpperOnly) {
      for (a = 0; a < NGauge; ++a)
        for (b = 0; b < NGauge; ++b)
          MT(a,b) += DebyeGauge[a][b] * std::pow(T, 2);
    } else {
      for (a = 0; a < NGauge; ++a)
        for (b = a; b < NGauge; ++b)
          MT(a,b) += DebyeGauge[a][b] * std::pow(T, 2);
      SymmetrizeLower(MT);
    }
  }
  return {needZero ? std::move(M0) : MatrixXd{}, std::move(MT), needZero};
}
```

`BSMPT_SAFE_GAUGE_UPPER_ONLY` 目前默认关闭。若该环境变量关闭，当前代码
确实分别计算上下三角，不能为了“矩阵对称”而改变为上三角，否则虽然数学
结果相同，浮点求和顺序和极端点行为都会改变。若打开该变量，pair builder
也必须保持当前的上三角后复制策略。

## `diff` 分支如何接入

`V1Loop` 开头根据分支决定是否需要零温矩阵：

```cpp
const bool needHiggsZero = !C_UseParwani;
const bool needGaugeZero = true;
auto H = BuildHiggsMassMatrixPair(v, Temp, needHiggsZero);
auto G = BuildGaugeMassMatrixPair(v, Temp, needGaugeZero);
```

### `diff == 0`

在 pair 的两个矩阵上分别执行现有的
`SelfAdjointEigenSolver<MatrixXd>(M, EigenvaluesOnly)`，并保留当前
`abs(eigenvalue) < ZeroMass` 时置零的逻辑。

- Parwani：使用 `H.thermal`、`G.thermal`，另外使用 `G.zero`；不使用
  `H.zero`。
- 非 Parwani：使用 `H.zero` 计算 CW 与 Daisy 的 Higgs 项，使用
  `G.zero` 计算 gauge 项；`H.thermal`/`G.thermal` 用于热质量项。

### `diff > 0`

当前的导数矩阵 `dM/dphi_i` 不含 Debye 温度项，因此同一 `dM` 可复用：

```cpp
MatrixXd dH = BuildHiggsDerivativeMatrix(v, diff);
MatrixXd dG = BuildGaugeDerivativeMatrix(v, diff);
```

分别对 `(H.thermal,dH)` 与 `(H.zero,dH)`、以及 `(G.thermal,dG)` 与
`(G.zero,dG)` 调用原有 `FirstDerivativeOfEigenvalues`。这仍然是两次原有
导数求解；仅减少矩阵组装。Parwani 分支只需 H(T) 的一套导数和 G 的两套。

### `diff == -1`

温度导数矩阵必须分别保留：

```cpp
dH_T = 2 * Temp * DebyeHiggs;
dG_T = 2 * Temp * DebyeGauge;
dH_0 = ZeroMatrix;
dG_0 = ZeroMatrix;
```

对零温质量矩阵调用原有 `FirstDerivativeOfEigenvalues(M0, ZeroMatrix)`，不
能把 `dM_T` 误套到 `M0`。这尤其重要，因为 `V1Loop` 的 `diff==-1` 会把
温度导数传给 `boson/fermion`，改变物理输出。

## 浮点顺序与准确性约束

为了让已有 TSV 对照尽可能保持 bitwise 一致，实施时应遵守：

1. `L2 -> L3(k) -> L4(k,l)` 的循环顺序不变；不要改成矩阵乘法、BLAS
   rank update 或先按幂次重新分组。
2. `M0` 先完成，再以 `MT = M0` 复制；Debye 项仍按原矩阵索引和
   `+= Debye * std::pow(T,2)` 顺序加入。
3. 保持 `MatrixXd` 与原 `SelfAdjointEigenSolver<MatrixXd>`；不切换固定维度
   类型，不切换 eigensolver，不缓存 eigenvectors。
4. 保留当前 `ZeroMass` 阈值和 eigenvalue 顺序；不要为接近简并的谱做排序
   或重标号。
5. 质量矩阵的 pair 只复用组装结果；每个 `M0` 与 `MT` 仍各自运行原有
   eigensolver。这样不会把“矩阵复用”误变成“本征值复用”。

严格 bitwise 一致仍需用 `compare_outputs.py --rtol 0 --atol 0` 在至少一个
`mH<125` 和一个 `mH>125` 点上检查；若只达到数值容差而非 bitwise，应记录
差异，不能默认为可接受。

## 缓存生命周期与多输入行安全性

推荐缓存生命周期为**单次 `V1Loop` 调用**：

```text
进入 V1Loop(v,T,diff)
  构造局部 H/G pair
  完成本次势能的所有 H/G 质量和导数使用
退出 V1Loop，释放 pair
```

不建议把 `(v,T)` 矩阵放入 `Class_Potential_Origin` 的 mutable 成员缓存：

- CalcGW 的 bounce 栅格会连续传入不同 `v` 和不同 `T`；跨调用命中率低，
  而 8x8/4x4 矩阵缓存会产生大量比较和内存流量。
- `set_gen` 是派生类 virtual 函数，外部可以直接调用它；仅在基类
  `resetbools()` 中清缓存不能覆盖所有模型参数变化路径。
- 不同输入行复用同一个模型对象时，`Curvature_*`、`Debye*`、scale、
  `NHiggs/NGauge` 等都可能变化；旧 cache 若未完整失效会污染下一行。
- CalcGW 可能由外层多进程/线程使用；共享 mutable cache 还会制造数据竞争。

局部 pair 不依赖行号、参数存储或全局状态，因此天然满足：每一行新建/设置
的模型参数只影响本次 `VEff`；下一行不会看见上一行矩阵。如果未来确实要
做跨 `V1Loop` cache，必须至少把完整 `(v,Temp,diff,scale,parStored,
parCTStored,SetCurvatureDone,Debye*)` 版本号纳入 key，并在所有模型的
`set_gen`/`set_All`/`resetScale` 路径失效；这不再是最小改动，也不建议作为
第一版优化。

## 验证方案

先实现仅 `diff==0`、非 Parwani 的 H/G pair，使用两份独立二进制：

1. 逐个比较 `HiggsMassMatrix(v,T)` 原矩阵与 pair 的 `thermal`，以及
   `HiggsMassMatrix(v,0)` 与 pair 的 `zero`；使用 element-wise max abs/rel
   误差和 bitwise 检查。
2. 比较 `HiggsMassesSquared` 生成的质量向量与 pair + 原 solver 的向量。
3. 比较完整 CalcGW TSV，`--rtol 0 --atol 0`；分别覆盖 `mH<125`、
   `mH>125`，并包含极小 SNR 与较大 SNR 点。
4. 再打开 Parwani 编译/运行测试（即使生产配置当前为非 Parwani），确认
   没有错误访问未构造的 `H.zero`。
5. 最后再接入 `diff>0` 和 `diff==-1`，分别测试场梯度和温度导数。

性能计时要单独报告矩阵组装计数和总 `V1Loop` 时间；由于 eigensolver 不变，
若总耗时没有明显变化，说明瓶颈在本征值求解、热函数或 bounce，而不是继续
把矩阵循环微优化。

## 风险边界

本设计不包含以下改动：

- 不改变任何矩阵维度或块对角化假设；
- 不替换 Eigen eigensolver；
- 不改变热函数积分/插值迭代；
- 不降低 bounce 栅格、积分精度或收敛阈值；
- 不修改 thdmtools、原项目或系统连接。

因此它适合作为下一轮 CalcGW 精确模式的低风险代码级实验，但实际收益必
须以运行计数和 mH 两个区域的逐字段对照决定，而不是仅由理论 FLOP 数估计。

## 后续只读审查：同一调用中的其它可复用项

以下结论来自继续检查 `V1Loop`、`VEff`、GW 后处理和热函数实现；本节同样
没有修改源码。

### V1Loop 中已经只算一次、无需额外 cache 的量

在 `ClassPotentialOrigin.cpp:3112-3118`，每次 `V1Loop` 已经只调用一次
`QuarkMassesSquared(v,diff)` 和 `LeptonMassesSquared(v,diff)`。它们与温度无关，
后面的各个分支直接复用向量；没有发现再次组装夸克/轻子质量矩阵的路径。
非 Parwani 分支中的夸克颜色简并也已经先累计 `AddContQuark` 再乘
`NColour`（约 `3164-3169`、`3240-3247`、`3350-3354`），不应重复优化。

`HiggsMassesVec`、`GaugeMassesVec` 及其零温向量在各分支中也已经以向量形式
重复使用。真正的重复仅是生成这些向量之前的 H/G 矩阵与 eigensolver；这正是
本文前半部分的 pair 方案所针对的目标。

### `Temp == 0` 时可安全复用整个质量向量

在 `V1Loop:3114-3116`，当 `Temp == 0` 时：

```text
HiggsMassesSquared(v, Temp, diff) == HiggsMassesSquared(v, 0, diff)
GaugeMassesSquared(v, Temp, diff) == GaugeMassesSquared(v, 0, diff)
```

这对 `diff==0`、`diff>0` 和 `diff==-1` 都成立：温度矩阵的 Debye 项为零，
温度导数矩阵 `2*Debye*Temp` 也为零。可以在 `Temp == 0` 的极少数调用中把
零温向量复制给热向量，避免第二次 eigensolver；这不是当前 CalcGW 主要路径
（通常 `T>0`），应作为独立的边界特例测试。为了保持现有代码的调用顺序，
第一版仍建议只复用矩阵而不改变 eigensolver 数量；若启用该特例，必须以
bitwise TSV 对照确认 Eigen 对同一矩阵的第二次求解没有产生可观察差异。

### `Getkappa_col` 存在跨函数的同点重复

`src/gravitational_waves/gw.cpp:708-709` 对 false/true vacuum 各调用一次
`VEff`。`VEff` 在 `ClassPotentialOrigin.cpp:3091-3095` 又调用完整 `V1Loop`，
其中已生成一遍 H/G/Q/L 质量谱。紧接着 `gw.cpp:721-737` 又对两个完全相同
的真空和 `Tstar` 调用 H/Q/L/G 质量函数。因此在最终 GW 后处理中，质量矩阵
与 eigensolver 会再做一遍；这不是跨输入点 cache，而是同一函数内可以消除的
结构性重复。

低风险的长期方案是增加一个局部返回结构（不放入模型成员）：

```cpp
struct VEffWithSpectrum {
  double oneLoopValue;
  std::vector<double> higgs, gauge, quark, lepton;
};

// 对一个 (v,T,diff) 只组装/求解一次质量谱，并返回 V1Loop 和 spectrum
VEffWithSpectrum EvaluateVEffAndSpectrum(
    const std::vector<double>& v, double T, int diff = 0) const;
```

`VEff` 可以只返回该结构中的标量；`Getkappa_col` 则使用 false/true 各一个
局部结果，同时取得 `dV` 所需的 NLO 势和压力所需的质量谱。结构只在
`Getkappa_col` 调用栈内存在，输入行间、线程间不会共享。实现上应先确认
`V1Loop` 的质量向量已经按当前分支计算完，再将现有的 `res` 与向量一并返回；
不应在外层引入全局 map 或按 `v` 的哈希缓存。

这是比 H/G pair 更大的一项 API 改动，适合作为单独实验；验证时必须比较
`dV`、`P_LO`、`P_NLO`、`kappa_col` 和完整 GW TSV，而不只比较质量向量。

### `V1Loop` 的温度常数可局部预计算，但不保证值得做

`ClassPotentialOrigin.cpp:91-124` 的 `boson` 和 `:128-160` 的 `fermion` 在
每一个粒子项中重复计算 `pow(Temp,2)`、`pow(Temp,3)`、`pow(Temp,4)`、
`pow(M_PI,2)` 和质量比。这些值在一次 `V1Loop` 内相同；可以引入一个内部
`ThermalFactors{T2,T3,T4,inv2Pi2}`，或只在 `V1Loop` 内传入预计算因子。

但是这会改变乘法/舍入顺序，不能宣称 bitwise 等价；而且总共约几十个粒子
项，耗时预期远小于 H/G eigensolver。建议仅在 profiler 证明热函数占明显比例
时实验，并使用数值容差和完整 TSV 对照。

### 热函数没有运行时积分迭代可删除

`src/ThermalFunctions/ThermalFunctions.cpp:74-199` 的正区间实现是固定阶数
低/高区间展开（低阶循环最多 `l=2..5`，高阶最多 `l=0..3`）；
`:216-230` 只做区间分支；负玻色子在 `:235-247` 直接查预构造 cubic spline。
没有每个质量点调用 GSL 自适应积分或迭代求根，因此不能通过“减少积分迭代
次数”保持准确地加速这里。系数计算器在构造期间预计算并且文档标明为线程安全，
不需要再加共享 cache。

### `Getkappa_col` 的 spectrum 合并需要保留物理分支差异

`Getkappa_col` 的 `P_LO` 使用 H/G/L/Q 全部质量平方，而 `P_NLO` 只使用有序
的 gauge 质量（`gw.cpp:741-782`）。合并 API 时不能只返回 H/G，也不能以
`VEff` 内部使用的 Daisy/Parwani 质量替代外层明确请求的
`HiggsMassesSquared(v,T)` 等谱。当前生产配置是非 Parwani，但设计必须让
`C_UseParwani` 的两套语义分别保留；否则 `kappa_col` 会出现静默物理改变。

## 推荐实施顺序（仅作为主 agent 后续实验输入）

1. 先实现本文前半部分的局部 H/G pair，保持原 solver 数量；这项改动最小，
   可直接用两个质量区域的已有 CalcGW TSV 验证。
2. 单独计数 `Temp==0` 特例的命中率；若 CalcGW 没有命中，暂不加入生产路径。
3. 用 profiler 测量 `Getkappa_col` 中 VEff 后再次求谱的耗时，再决定是否实现
   `VEffWithSpectrum` 合并；该项不能与 H/G pair 同时首次引入，便于定位差异。
4. 最后才考虑热函数常数预计算；它属于微优化，优先级低于质量矩阵复用。

## 对已实现 Higgs pair 的 `diff>0` 扩展审查

当前 opt-in 实现位于 `ClassPotentialOrigin.cpp` 的 `V1Loop`（约
`3210-3270`，实际行号随其它隔离实验变化），条件严格限定为 R2HDM、
`diff==0`、非 Parwani。继续检查 `HiggsMassesSquared`（约
`2478-2560`）后，`diff>0` 的扩展在数学上也有一个安全复用点：

```text
M_H(v,T)       = M_H,field(v) + T^2 DebyeHiggs
M_H(v,0)       = M_H,field(v)
dM_H/dphi_i   = L3[...,i] + sum_k L4[...,i,k] v[k]
```

`HiggsMassesSquared(v,Temp,diff)` 先组装温度质量矩阵，再独立组装上述
`dM_H/dphi_i`；随后 `HiggsMassesSquared(v,0,diff)` 又重新组装同一个质量
导数矩阵。对同一 `V1Loop` 调用，`dM_H/dphi_i` 不含 Debye 项，所以可以只
组装一次，并将同一个数值矩阵传给两个独立的
`FirstDerivativeOfEigenvalues`。

建议采用独立开关进行实验，例如：

```cpp
const bool use_r2hdm_higgs_pair_diff =
    env("BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF") &&
    Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 && nVEV == 4 &&
    !C_UseParwani && diff > 0 && diff <= static_cast<int>(NHiggs);
```

暂不把它合并到现有 `BSMPT_USE_R2HDM_HIGGS_PAIR`，以便分别测量 `diff==0`
和梯度路径的差异。

### 严格保持两次导数 eigensolver

不能因为 `dM` 相同而复用 eigenvectors 或导数结果：

- `M_H(v,T)` 与 `M_H(v,0)` 一般不同，Debye 项会改变 eigenvectors、简并分组
  和 `FirstDerivativeOfEigenvalues` 的 mapping。
- 必须分别对 `M_T`、`M_0` 调用原始
  `FirstDerivativeOfEigenvalues(MassCast, DiffCast)`，维持两次
  `SelfAdjointEigenSolver<MatrixXcd>`。
- 允许复用的是 `Diff` 的多项式组装；为最严格的 bitwise 对照，可以把共享
  `MatrixXd Diff` 分别复制成两个 `MatrixXcd DiffCast`，不共享可变 solver 状态。

建议伪码：

```cpp
MatrixXd M0 = HiggsMassMatrix(v, 0);       // 原有循环顺序
MatrixXd MT = M0;
AddDebyeToUpperTriangleAndMirror(MT, Temp);

MatrixXd dM(NHiggs, NHiggs);
for (i = 0; i < NHiggs; ++i)
  for (j = 0; j < NHiggs; ++j) {
    dM(i,j) = Curvature_Higgs_L3[i][j][diff - 1];
    for (k = 0; k < NHiggs; ++k)
      dM(i,j) += Curvature_Higgs_L4[i][j][diff - 1][k] * v[k];
  }

MatrixXcd dM_T = dM;                       // 独立输入对象
MatrixXcd dM_0 = dM;
MatrixXcd mT = MT;
MatrixXcd m0 = M0;
auto thermal = FirstDerivativeOfEigenvalues(mT, dM_T);
auto zero    = FirstDerivativeOfEigenvalues(m0, dM_0);
```

`dM_T` 与 `dM_0` 的复制不是性能重点，但可排除未来 solver 修改导致的
aliasing 风险；两个 `FirstDerivativeOfEigenvalues` 本身仍完全独立。

### 与当前 `V1Loop` 分支的接入边界

当前非 Parwani `diff>0` 分支在初始质量谱之后，于约 `3227` 再计算
`HiggsMassesSquared(v,0,diff)`；这是唯一需要 H(0) 导数谱的分支。Parwani
`diff>0` 只消费 H(T) 导数谱，不应为了 pair 实验额外构造 H(0)。因此扩展
开关应满足：

```text
非 Parwani + diff>0: pair(MT,dM) + pair(M0,dM)
Parwani    + diff>0: 原 H(T) 路径
diff==0              : 现有 diff==0 pair 或原路径
diff==-1             : 原路径（dM/dT 在 T 与 0 不相同）
```

为防止未初始化向量，进入非 Parwani分支时应以 `use_pair_diff` 判断是否跳过
原来的 H(0) 调用；并继续保留 `HiggsMassesVec.size()` 检查。

### 可安全复用与不可复用的工作清单

| 工作 | `diff>0` 是否可复用 | 原因 |
|---|---|---|
| H 场依赖多项式矩阵组装 | 是 | T 只增加 Debye 项 |
| H `dM/dphi_i` 组装 | 是 | Debye 与场无关 |
| H(T)/H(0) eigenvalues | 否 | 两矩阵不同，必须两次 solver |
| eigenvectors / degeneracy mapping | 否 | 由各自矩阵决定 |
| `MatrixXd -> MatrixXcd` 导数转换 | 可复制共享数值 | 但建议保留两份输入对象以保证隔离 |
| `boson(..., diff)` 值 | 否 | 质量平方不同，且 CW/热函数逐项输入不同 |

### 线程安全与输入行安全

扩展必须继续使用 `V1Loop` 栈上的 `M0`、`MT`、`dM`；不能把导数矩阵放入
模型成员。这样：

- 每个 `V1Loop` 调用有自己的矩阵和两个 solver 输入；
- 两个 CalcGW 进程或线程不会共享矩阵；
- 下一输入行通过新的 `V1Loop` 局部对象自然清空状态；
- 即使外部直接改变模型参数，旧调用栈之外也不存在 cache。

静态环境变量只用于读取开关，不保存质量矩阵；其初始化由 C++ 保证线程安全。
应在每个进程启动前设置环境变量，不支持运行中动态切换。

### 严格验证顺序

扩展实现后应按以下顺序比较：

1. 对固定的 `v,T,diff=1..8`，比较原 `HiggsMassesSquared(v,T,diff)` 与 pair
   的 thermal 导数向量，及原 H(0) 与 pair 的 zero 向量。
2. 至少覆盖一个 `mH<125` 和一个 `mH>125` 输入点；先单独运行梯度探针，再
   运行完整 CalcGW。
3. 用 `compare_outputs.py --rtol 0 --atol 0` 比较非 runtime 字段；若
   bitwise 不一致，再报告 max abs/relative，而不是直接开启生产扫描。
4. 开关关闭时必须确认二进制输出回到当前基线；Parwani、`diff==0` 和
   `diff==-1` 的路径不应受该扩展影响。

目前不建议直接扩大现有开关的条件；先以独立开关、单点调用计数和两质量区域
对照验证，确认导数复用没有在近简并 eigenvalue 或极小 SNR 点引入状态差异。
