# 非严格 CalcGW 专项优化与执行方案

更新日期：2026-09-05。本文件只描述项目主目标。概念、实现、验证和停止条件均在此
维护；不得再建立平行的“当前非严格路线图”。

实施任一候选前必须查`RESEARCH_INDEX_ZH.md`的危险证据和禁止重复表，并在
`RESEARCH_JOURNAL_ZH.md`登记新轮次；后续AI的工作纪律见`AGENTS.md`。

## 1. 当前问题

现有首遍 central2 + adaptive64 + raster500 + thermal-fast 在旧42行矩阵中由
2173.759秒降到约908.155秒，首遍约快58.2%，但默认guard只接受6行、回退36行。
大量点实际支付`approx + exact`，所以继续压缩首遍几个百分点不能解决主问题。

现有guard依赖输出状态、极弱SNR和少数风险锚点。它能捕获已知假阳性，但无法在运行
严格版前估计一个状态看似正常的近似结果是否已有10%以上SNR偏差；classified E1先导
中的54.95%反例证明了这个缺口。

## 2. 目标架构

每点按以下级联执行：

```text
BSMPT input
  -> R0 覆盖/OOD检查
  -> L1 当前低精度近似
  -> C1 廉价内部证书
       -> 高置信安全：接受L1
       -> 可修复不确定：L2局部加密/中精度
                            -> C2通过：接受L2
                            -> C2失败：exact-fast
       -> 明确危险或OOD：exact-fast
```

R0只能使用BSMPT输入和纯BSMPT验证覆盖；C1/C2只能使用当前计算产生的状态、残差和
双分辨率差异。禁止读取thdmTools结果、历史严格SNR或文件类别。

## 3. N1：只记录内部诊断，不改变算法

目的：找出能够在严格重算前预测状态翻转或SNR超限的量。首先在近似overlay增加实验
输出或独立日志，不改变任何返回值：

- action温度节点、相邻斜率和插值midpoint误差；
- bounce shooting的undershoot/overshoot括区间宽度、迭代数和终端残差；
- ODE末端边界残差、profile单调性、virial/action一致性；
- path deformation法向力、迭代收敛率和路径跳跃；
- raster上最大局部曲率、符号变化、最小质量平方绝对值和thermal分段距离；
- `Tn/Tp/T*`定位所跨节点数，`beta/H`使用区间宽度和局部拟合条件数；
- coarse结果距离SNR=10/100等业务阈值的相对余量。

实施：

1. 用4个正常点确认诊断不改变非runtime输出；
2. 用已知A/B反例、B扰动云和54.95%反例验证诊断是否异常；
3. 在旧42行及E1辅助点上保存诊断与真实误差的联合表；
4. 只保留计算增量低于L1耗时2%、且能明显区分反例的指标。

N1本身不扩大接受面，也不宣称安全区。若已知54.95%反例与安全点在所有候选指标上
不可区分，则单纯输出证书路线失败，进入N2主动双分辨率而不是训练黑箱猜测。

### N1具体源码落点

所有编辑只发生在实验副本 `bsmpt_speed_lab/upstream/`，并以环境变量默认关闭：

- `src/bounce_solution/calcgw_profiler.cpp`及对应header：扩展独立诊断计数器和JSON/TSV
  sidecar输出；不开启时不得增加热路径原子操作。
- `src/bounce_solution/action_calculation.cpp::RasterizedVdl`：记录500格局部曲率、符号
  变化和候选中点索引；N1不额外调用VEff。
- 同文件`Solve1DBounce`及threshold搜索：记录括区间、迭代数、终端误差、最终profile
  边界量和adaptive threshold选择，不改变搜索顺序。
- `src/bounce_solution/bounce_solution.cpp::CalculateActionAt`与GW高/低action扫描：记录
  `SolutionList`温度/action节点、相邻斜率、插值区间和后续`beta/H`使用区间。
- `run_calcgw_approx_*`：新增仅研究用诊断开关和sidecar路径；严格wrapper不传递该变量。

N1a首先只加入已有量的记录；需要额外VEff/action求值的内容全部属于N2，不得混入N1
的零扰动和成本基线。

## 4. N2：低成本双分辨率证书

不直接再跑完整raster1000，而对L1结果做嵌套、可复用的局部探针：

### N2a action温度独立midpoint

在决定`Tn/Tp/T*`和`beta/H`的关键插值区间额外计算未参与拟合的midpoint。比较真实
action与预测值、根位置移动和局部导数。超过阈值时只细分该区间，不重算全温区。
实现位于`bounce_solution.cpp`：在原`SolutionList`已排序后选取真正决定根和导数的区间，
调用现有`CalculateActionAt(midpoint)`，随后重新使用原有spline/温度求解路径。探针结果
必须标记，禁止反过来作为同一轮证书的预测训练节点。

### N2b raster局部加密

从500格结果选择高曲率、符号变化、零点和质量阈值附近区间，插入独立中点；复用已有
端点。比较`dV/dl`预测误差、action变化和射击终点。异常才升级相关区间或exact。
实现位于`action_calculation.cpp::RasterizedVdl`：保留原500格数组，只对预先选择的区间
计算中点真实`Calc_dVdl`；若升级则构造新的有序数组再调用同一`spline::set_points`。

### N2c bounce复核

对最终候选profile使用更严ODE步长或独立shooting bracket作一次廉价复核，比较action、
末端边界和profile形状。不得用相同离散数据同时拟合和认证。
实现继续复用`Solve1DBounce`，但使用独立、更严配置创建复核对象；不得原地污染L1的
profile、spline hint或SolutionList。C++对象生命周期先限定在单参数点，禁止跨行缓存。

阈值只可用E2训练部分选定；必须同时设置“安全、未知、危险”三段，禁止强迫每点二元
分类。证书成本目标不超过exact-fast的10%，否则难以在当前回退率下获得净收益。

## 5. N3：有证书的算法加速

只有N1/N2能捕获全部已知反例后才实施：

1. profile/shooting continuation：相邻温度复用profile、括区间和尺度；旧括区间必须
   重新验证包根，残差不降即冷启动；
2. safeguarded Brent/Illinois：始终保留两侧，候选出界或残差不连续即回二分；
3. adaptive action/raster：由N2同一独立探针决定局部细分，不能建立证书则恢复完整；
4. SpectrumJet方向导数：简并、过零和thermal分段附近以原差分抽查并强制回退。

每项单独比较；不得把continuation、求根和adaptive首次捆绑，以便定位错误来源。

## 6. N4：覆盖与风险路由

router只负责减少明显不适合L1的浪费，不负责替代数值证书：

- 输入超出E2/E3覆盖凸包、Yukawa/模式未覆盖或离最近验证点过远：直接exact；
- 局部邻域中出现任一unsafe点：扩大缓冲区或direct exact；
- 仅当局部E2支持且C1/C2通过，才接受近似；
- 模型优先采用可审计的标准化kNN距离/局部椭球和单调阈值，不先使用复杂分类器；
- archive category、source身份和所有`th_*`字段禁止成为生产特征。

E1只用于提出候选和找反例，不能扩大R0覆盖。router阈值冻结后才生成E3。

## 7. 样本设计

### E1辅助阶段

复用已有严格/近似配对和逐行认证的classified严格结果。分层覆盖BSMPT状态、SNR量级、
参数网格、source和机器。其通过率只写作“E1样本内”，不能报告总体概率。

### E2定向阶段

每个候选安全锚点和反例至少生成轴向`±δ`、随机方向扰动及安全/不安全连线二分点。
覆盖单/多步切换、NLO、共存相、bounce、弱信号、SNR阈值和质量简并。完全不运行
thdmTools，每点严格/近似双跑。

### E3确认阶段

冻结输入域、生成分布、分层比例、随机种子、router和证书阈值。首轮至少300个有效
独立点，按参数空间块留出；若目标是约95%置信下错误率低于1%，需零失败约299点。
任何失败都要报告实际区间并收缩规则，不能删点后重新计算成功率。

## 8. 正确性标签

`safe`要求：

- 所有共同status、活动相变块数、有限值掩码、`transition_history`一致；
- 每个有效总SNR和col/sw/turb分量相对误差≤10%；
- 关键温度、alpha、beta/H通过预先冻结的绝对/相对容差；
- 输出来自真正近似路径，不是exact fallback。

严格SNR接近零、不可追溯配置和覆盖外点为`unknown -> exact`。任一状态翻转、假/漏信号
或SNR超限为`unsafe`，优先级高于平均误差和速度。

## 9. 性能核算与门槛

逐层记录：R0耗时、L1耗时、C1耗时、L2触发率/耗时、exact回退率/耗时和峰值RSS。

```text
T_total = T_R0 + T_L1 + T_C1
        + P(L2) * (T_L2 + T_C2)
        + P(exact) * T_exact
gain = 1 - T_total / T_exact
```

晋级门槛：

- 已知反例和E2测试零漏放；
- 接受点状态/history零错误、SNR分量≤10%；
- 7%–10%点以当前严格版重复确认；
- E2组总时间至少快2%，E3目标快20%；
- 若高回退使总时间不优于direct exact，立即淘汰或收缩L1适用域。

## 10. 固定实施批次

1. **N1a**：为当前approx overlay加入只读计数/诊断，4点零扰动验证。
2. **N1b**：A/B、B扰动云、54.95%反例，检查零漏放潜力。
3. **N1c**：旧42行建立诊断—误差表，估计成本和过度回退可减少量。
4. **N2a**：只实现action关键区间midpoint证书；失败不进入其它N2。
5. **N2b/N2c**：分别测试raster局部加密和bounce独立复核。
6. **E2-1**：围绕所有已知unsafe及候选safe生成纯BSMPT扰动。
7. **N3**：按continuation→shooting→adaptive→SpectrumJet顺序单变量实现。
8. **N4**：冻结覆盖/OOD规则，先只减少必然回退的近似首遍。
9. **E3**：规则冻结后的独立确认与端到端统计。

每批完成立即更新`RESEARCH_JOURNAL_ZH.md`并提交。最低阶梯失败即停止该候选，不通过
不断追加规则掩盖反例。
