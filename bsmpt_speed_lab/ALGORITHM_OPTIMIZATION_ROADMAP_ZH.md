# CalcGW 算法级优化路线图

更新日期：2026-09-05。本文取代“继续依靠PGO/链接参数获得主要收益”的思路。
所有实现仍只能位于 `bsmpt_speed_lab/`；主项目、冻结严格分支和严格二进制不修改。

历史 classified 数据的证据等级、选择偏差和禁止用途由
`BSMPT_OPTIMIZATION_AND_BIASED_DATA_PROTOCOL_ZH.md` 约束。该数据只能作为E1辅助
回归，不能定义完整参数空间、安全路由或总体错误率；任何安全区域必须经不经过
thdmTools的纯BSMPT E2/E3样本确认。

## 1. 已确认的复杂度账本

代表性严格点约35秒中，phase tracing约8.3秒，bounce约26.9秒；其中24次
`Solve1DBounce` 合计约21.5秒。一次4维 numerical gradient需要16次 VEff，
Hessian需要33次，combined gradient+Hessian需要41次。每次 action 的1001点
raster约触发16016次 VEff，dense threshold约触发41041次 VEff。

现有代码已经用 `WarpPath` 复用相邻温度的路径，并在新温度选择最近的成功路径；
因此“加入路径warm start”不是新优化。尚未复用的是径向 bounce profile、射击括区间、
最佳 `l0/L` 和同一温度路径变形轮次间的一维解。

## 2. 严格算法路线 S1–S6

### S1：消除确定性重复求值

首个单变量实验是 `MinimumTracer::LocateMinimum`：循环条件和循环体在相同
`new_guess` 上重复调用 `df`。缓存该次gradient，不改变gradient内部公式或求和顺序。
随后逐项审计 BackwardsPropagation 收敛点Hessian和关闭日志下的诊断导数；历史已有
helper或实验开关的项目必须先确认是否已经验证，禁止重复声称新收益。

验收：除runtime外逐字段零差异；高SNR两次收益低于2%则不扩展。预计整体3%–15%，
但以实测为准。

### S2：dense threshold 确定性并行

1001个候选threshold点的MinSol、LinSol和Error彼此独立。固定输入snapshot后由最多4个
worker按索引计算，写入预分配槽位；主线程严格按原索引顺序选择MinError。禁止并行
reduction、共享spline hint和worker日志。

验收：先用2线程，Eigen/BLAS/OpenMP固定1线程；逐字段零差异并记录峰值RSS。超过
18GB、出现竞态/TSAN问题或总体收益低于5%即停止。预期该子阶段2–3.5倍。

### S3：RasterizedVdl 确定性并行

1001个 `l[i]`/`dVdl[i]` 独立计算，worker只写固定索引；主线程按顺序调用
`set_points()`。`Spline`只读，`RasterizeddVdl`、hint、status和Logger禁止worker写入。
S2和S3分别测试，不能首次捆绑。

验收同S2。预期raster子阶段2–4倍，整体可能1.2–2倍；实际占比低于10%则停止。

### S4：potential action map + ordered reduce

约2000个积分区间先并行计算各自term，主线程仍按原顺序累加，避免并行reduction改变
末位。预计整体5%–20%，但当前profile中action integral仅约0.66秒，优先级低于S2/S3。

### S5：path deformation内部并行

只并行同一轮内约500个path check采样或各knot force；迭代轮次保持串行，最大值和
收敛判定按原索引串行归并。预计子阶段1.5–3倍，但当前占比约1.8秒，后于S2/S3。

### S6：coexistence pair粗粒度并行

phase数据冻结后，不同coexistence pair可各自构造局部BounceSolution，完成后按原pair
编号归并。只对至少两个昂贵pair启用，最多2个worker，且与action内部线程池互斥。
单pair点无收益；多pair预计1.5–2倍。

## 3. Guarded算法路线 G1–G5

### G1：SpectrumJet / Hessian-vector

为完整R2HDM质量矩阵构造 `SpectrumJet`，一次谱分解同时产生V1Loop、四维gradient和
路径方向Hessian。使用Hermitian矩阵导数与divided difference处理简并子空间，保留
完整12×12 quark和9×9 lepton矩阵，不重复已失败的6×6降维。

第一步只替换路径方向二阶导数，原16次VEff gradient保持不变，单点约从41次VEff降到
17次；通过后才扩展完整jet。预计混合版整体1.3–1.8倍，完整jet可能2–4倍。

强制回退条件：本征值间隔过小、质量靠近零/CW/Daisy截断、thermal分段边界、非有限
矩阵，或独立3/5点方向差分校验失败。它不进入冻结严格入口。

### G2：径向profile与射击括区间continuation

相邻温度保存 `rho_sol/l_sol/dldrho_sol`、undershoot/overshoot括区间和最佳`l0/L`；
新温度先验证旧括区间仍包根，再尝试warm start，失败立即执行原冷启动。必须检查ODE
残差、末端边界、单调性、法向力和action virial。薄壁、Tc附近、势垒消失、多步跳跃
和action非单调点直接冷启动。预计action加速20%–60%。

### G3：有括区间保护的Brent/Illinois射击

保留undershoot/overshoot两侧，以连续终端残差做safeguarded求根；残差不连续、候选
出界或单调性破坏就退回二分。预计整体5%–25%，必须在G2独立完成后单变量测试。

### G4：带后验误差的adaptive raster/action温度节点

不再固定把1000降为500。端点/中点真实求值与插值预测比较，在曲率、符号变化、零点、
质量阈值附近递归细分；无法建立误差证书就完整回退1000点。action温度采样同时比较
线性与单调PCHIP，并使用独立midpoint hold-out。预计减少30%–60% action节点或
50%–85%平滑raster求值。

### G5：Anderson/L-BFGS路径变形与多场BVP

短期可对Bernstein系数残差使用带信赖域的Anderson/L-BFGS，法向力增大立即恢复旧路径；
长期可研究直接多场collocation/BVP continuation。前者预计整体10%–30%，后者在平滑
序列可能2–10倍，但风险和验证成本最高，最后实施。

## 4. Phase算法路线 P1–P4

- P1：`Phase::Get` 线性查找改为有序 `lower_bound`，严格但预计低于3%。
- P2：predictor-corrector phase tracing，corrector仍用原求解；残差不降、Hessian异常、
  minimum跳跃即回退，属于guarded。
- P3：把 `CalculateTc` 的100点诊断扫描与求根需求分离；可能改变lazy phase填充顺序，
  只能先作为guarded实验。
- P4：Mode2每轮low/high side并行，完成后串行更新gap；局部预计1.3–1.8倍。

## 5. 线程与内存纪律

- 外层最多2个CalcGW；测试内部并行时优先只运行1个CalcGW。
- action worker从2开始，最多4；只有RSS证据充分才测试8。
- pair并行与action并行不得同时启用。
- Eigen、OpenMP、BLAS均固定1线程；现有minimizer多线程同时关闭。
- worker禁止写Logger、spline、SolutionList、status和共享hint。
- 峰值RSS达到18GB立即停止，给WSL 24GB上限保留至少6GB。

## 6. 固定实施顺序

1. S1 LocateMinimum重复gradient。
2. 为S2/S3增加纯计数/计时，确认threshold与raster真实占比。
3. S2 dense threshold两线程确定性并行。
4. S3 raster两线程确定性并行。
5. S4/S5及多pair样本上的S6。
6. G1混合SpectrumJet方向Hessian。
7. G2 profile/bracket continuation，再单独测试G3。
8. G4 adaptive误差证书。
9. G5路径加速；多场BVP仅作长期独立原型。

每轮保存调用次数、worker数、RSS、回退率、残差证书、原始TSV、全部status、history、
温度、alpha、beta/H、各SNR分量和墙钟。严格路线要求除runtime外逐字段零差异；
guarded路线要求状态/history一致、SNR分量误差不超过10%，7%以上复验，边缘点回退。
