# BSMPT Speed Lab agent instructions

本文件适用于 `bsmpt_speed_lab/` 整个子树。进入本目录工作的AI必须先完整阅读本文件，
再按任务路由读取最少必要文档。不得把历史报告误当作当前计划。

## 1. 永久边界

1. 只允许修改 `bsmpt_speed_lab/`；原版BSMPT、项目主程序、冻结严格分支/标签和
   `install/bin/CalcGW`不得修改。
2. 项目主目标是非严格CalcGW的安全、低回退、端到端加速；严格优化只是次要目标。
3. 最终方案只依赖BSMPT/CalcGW，不调用thdmTools，不用任何`th_*`字段或classified
   类别作优化、router、证书或正确性特征。
4. MultiNest classified是有偏E1辅助样本，不能覆盖完整参数空间；E1通过不能证明
   安全区，E1失败可淘汰候选。安全边界必须由纯BSMPT E2/E3建立。
5. 外层最多两个CalcGW；内部并行实验时只运行一个；构建最多`-j2`；RSS到18GB停止。
6. 每轮只改变一个独立变量。先做最小阶梯，失败立即停止；不得通过不断添加特例隐藏
   反例。成功和失败都更新日志并保存产物。
7. 非严格最高级标准是定性双向一致：`strict SNR>0`当且仅当`approx SNR>0`；严格
   fail/non-positive当且仅当近似fail/non-positive。假阳性或假阴性直接淘汰候选。
   SNR幅值10%只在双方positive后作为次级警戒，不能取代定性判断。
8. 旧guard 6/42接受、36/42回退不可作为目标方案。E3目标回退率≤20%、理想≤10%；
   新证书必须能够接受经验证稳定的失败结果，不能把所有failure无条件送exact。

## 2. 文档权威性

发生冲突时按以下优先级：

1. 本文件：agent行为、目录和任务路由；
2. `PROJECT_MASTER_PLAN_ZH.md`：项目目标、证据等级、成功定义和当前下一步；
3. `NON_STRICT_OPTIMIZATION_PLAN_ZH.md`：主目标的唯一当前执行方案；
4. `RESEARCH_INDEX_ZH.md`：已有研究和产物索引；
5. `RESEARCH_JOURNAL_ZH.md`：逐轮事实与commit；
6. `APPROX_SAFE_REPORT_ZH.md`、`OPTIMIZATION_REPORT_ZH.md`：历史技术事实；
7. `SNAPSHOT.md`：冻结strict快照。

README只导航。历史报告中的“推荐工作流程”“后续方向”如果与总纲或非严格专项冲突，
视为过时；不得据此重开已淘汰路线。

## 3. 任务路由

| 用户任务 | 必须阅读 | 下一步依据 |
|---|---|---|
| “继续/自动优化/接着做” | 总纲、非严格专项、日志最后一轮 | 执行专项“固定实施批次”中第一个未完成项 |
| 设计或实现非严格优化 | 非严格专项、研究索引、非严格历史报告相关节 | 先N1诊断，再N2证书；不得跳到静态区域分类 |
| 判断安全区域/router | 总纲E1/E2/E3、非严格专项N4/样本设计 | E1只提候选；E2建边界；规则冻结后E3确认 |
| 使用classified/MultiNest数据 | 总纲数据证据、研究索引数据部分、A22/A23日志 | 只认证真实CalcGW行，不用`th_*`；结果标E1 |
| 验证正确性/统计成功率 | 非严格专项正确性与性能、研究索引验证工具 | 分开报告E1/E2/E3和首遍/端到端结果 |
| 查询当前速度、回退率、错误 | 总纲当前状态、非严格历史报告、日志状态表 | 明确样本集和口径，禁止把42行推广到全空间 |
| 严格快速版优化 | 总纲次要目标、严格报告、快照 | 只考虑S1→S2→S3，逐字段零差异 |
| 复现某轮或寻找产物 | 研究索引、日志对应Ax、Git commit | 使用保存的TSV/wrapper，不凭摘要重造样本 |
| 构建/崩溃/double-free | 总纲资源规则、快照、严格报告安全注意 | 保持Eigen vectorization关闭，不覆盖正式安装 |
| 整理文档 | 本文件、README、总纲 | 当前计划只写总纲/非严格专项，历史事实只写日志/报告 |

## 4. 当前唯一下一步

默认任务是非严格专项N1a，而不是继续降低固定raster或训练参数空间分类器：

1. 先扩展比较器，分别输出严格、裸近似和guard最终的二类结果、假阳性、假阴性、
   最早失败阶段和回退原因；
2. 再在实验副本的`calcgw_profiler`、`RasterizedVdl`、`Solve1DBounce`、
   `CalculateActionAt`增加默认关闭、只记录已有量的诊断；
3. 用4个正常点证明诊断不改变任何非runtime输出且开销低于2%；
4. 用A/B及B扰动云检查所有定性翻转；54.95%点只作次级幅值检查；
5. 完成N1前，不实现N2额外求值、N3算法改变或N4上线router。

若用户明确改变主目标，先更新总纲和本文件，再执行；不得让新任务只存在于对话中。

## 5. 每轮工作模板

开始前：

- 查`RESEARCH_INDEX_ZH.md`确认是否已做过；
- 在日志声明候选、唯一变量、基线、数据等级和最低验证阶梯；
- 检查现有CalcGW进程和内存，确认输出只写本目录或`/tmp`。

完成后：

- 保存原始输入/输出、比较结果、runtime和RSS；
- 报告所有状态/history变化和每个SNR分量最大误差；
- 分开报告近似首遍、guard接受、严格回退和端到端时间；
- 更新日志的结论、产物、commit/push；
- 更新总纲“当前状态”及本文件“当前唯一下一步”，但只在阶段实际完成后更新。

## 6. 禁止重复

固定raster400/250/100、adaptive32、analytic-gradient、high-temperature fast powers、
boson coefficient pointer、exact-key thermal/VEff/mass cache、quark 12x12→6x6、
directional dV/dl、RK5 fixed/workspace、旧PGO复用、全程序当前PGO、容差键缓存、
quark dynamic noalias和无证书dynamic raster均已有失败或停止证据。除非出现新的机制、
明确不同的唯一变量和预先写明的重开条件，否则不得重复运行。
