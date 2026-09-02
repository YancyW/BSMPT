# CalcGW threshold 自适应搜索提案

本提案只针对 `src/bounce_solution/action_calculation.cpp` 的
`BounceActionInt::CalculateExactSolutionThreshold`，不涉及 thdmtools、数据筛选、
外层并行或原项目安装。补丁文件是 `adaptive_threshold_proposal.patch`，默认不启用。

## 为什么值得尝试

原实现把 `l0-lmin` 的对数区间固定划成约 1001 个点。每个点都会调用
`Calc_dVdl` 和 `Calc_d2Vdl2`，后者还包含数值 Hessian。因此 threshold 扫描会重复大量
势能、质量矩阵和本征值计算。提案先用 32 或 64 个粗点，再对多个局部候选区间加密。

## 使用方式

编译补丁后的隔离版时：

```text
BSMPT_ADAPTIVE_THRESHOLD=1 BSMPT_ADAPTIVE_THRESHOLD_GRID=64 CalcGW ...
```

`BSMPT_ADAPTIVE_THRESHOLD_GRID` 只接受 `32` 或 `64`；未设置时默认 `64`。环境变量
不存在、为空、以 `0` 开头，或网格值非法时，直接使用原来的约 1001 点扫描。

## 搜索算法与终止保证

1. 在原始指数区间端点（包含端点）取 33/65 个点。
2. 检查每个样本是否有限且可构成原有误差。粗网格局部极小值、两个端点以及按误差
   排序的前若干候选合并，最多保留 6 个候选。
3. 每个候选最多 6 层、每层 8 个子区间，仅在局部区间评价 7 个内部点。额外评价次数
   上限为 `6*6*7=252`，因此单次调用总预算有严格上界。
4. 若全部候选无效、出现非有限值/异常、指数区间无效，或粗网格出现极端相邻跳变，
   自动放弃自适应结果，执行未删减的原 1001 点扫描。
5. 原有 `MinError>1e-2` 时缩小 `FractionOfThePathExact` 的重试规则保留；提案用有限
   `while` 重现该递减过程，在 `FractionOfThePathExact<=1e-4` 时必然停止，不会因为
   fallback 递归嵌套。

## 准确性边界

这是保守原型，不应在未经对照时用于生产全量扫描。自适应搜索只能保证在已采样的
候选区间中寻找较小误差，不能数学上证明未采样区间没有更窄的极小值。因此建议对
`mH<125`、`mH>125` 分开，用 classified 中已有 CalcGW 结果逐点比较：SNR、相变温度、
作用量和所有非 runtime 字段均应在项目规定容差内一致；任何失败点都回退原扫描。

## 验证顺序

1. 先不设置环境变量运行基线，确认补丁没有改变默认行为。
2. 分别设置 `GRID=32`、`GRID=64`，对两个质量分支各选多个已计算点重复运行。
3. 用逐字段比较脚本检查输出，特别记录 threshold、action、SNR 及状态码；若发现
   质量分支或边界点偏差，停用自适应模式并保留原扫描。
4. 只有在准确性通过后，再统计 threshold 扫描调用次数和 CalcGW 总耗时。

补丁通过了针对当前隔离 upstream 源码的 `patch --dry-run` 检查；本轮没有修改
`bsmpt_speed_lab/upstream`、原 `install/bin/CalcGW` 或任何系统连接。
