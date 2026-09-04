# CalcGW A17 后续执行计划

更新日期：2026-09-04。本文把 A18–A20 固定成可执行、可停止、可反查的研究路线。
所有源码、构建、脚本和结果只能位于 `bsmpt_speed_lab/`；主项目、严格分支和严格
二进制禁止修改。编译最多 `-j2`，同时最多两个 CalcGW。

## 总体执行原则

1. Luna 只做只读审计：定位源码、核对历史、设计单变量实验和停止条件。
2. 主代理统一编辑、构建、运行、比较和 Git 操作，避免多个工作树互相覆盖。
3. 每个候选先做静态/构建资格检查，再做高 SNR 双配对；最低阶梯失败立即停止。
4. 严格候选要求除 runtime 外逐字段零差异。guarded 候选仍执行状态/history一致、
   每个 SNR 分量不超过10%、7%以上复验和危险边缘 exact fallback 的固定规则。
5. 单点收益必须重复两次，或五点合计至少2%；低于门槛不进入长矩阵。

## A18：选择性热点 PGO

### A18a：构建机制资格

- 将现有35个 `.gcda` 映射到源文件、对象和 CMake target。
- 查清 A17 中 `ClassPotentialOrigin_deprecated.cpp.gcda` 缺失是未执行、路径变化还是
  flags/重编译导致。
- 只给有完整画像的源文件增加 generate/use flags；未选择文件按 A8 control 编译。
- generate/use 必须保持对象路径、编译器、预处理宏和基础优化 flags 一致。
- 捕获完整构建日志；任何 selected source 的 missing、checksum mismatch、profile
  mismatch 或 partially trained warning 都判资格失败。

### A18b：热点集合

先测试最小集合，预计包括 bounce action、bounce solution、R2HDM有效势、热函数、
数值导数和 spline；实际集合以画像与调用热点证据交集为准。扩展集合只能在最小集合
已通过且有收益后加入，每轮只增加一组翻译单元。

Luna 双重审计后的 M0 固定为8个对象：`minimum_tracer.cpp`、
`ClassPotentialOrigin.cpp`、`bounce_solution.cpp`、`action_calculation.cpp`、
`transition_tracer.cpp`、`gw.cpp`、`const_velocity_spline.cpp`、`CalcGW.cpp`。
它们约占现有 time-profiler 累计热点的88%，且均有 current-system 画像。第一轮排除
`ClassPotentialOrigin_deprecated.cpp` 以及被多个 target 重复编译的
`MinimizeGSL/MinimizePlane/Minimizer/MinimizeNLOPT.cpp`。

实现采用实验室内 compiler-launcher，同时匹配源文件、输出对象路径和 target；不能
使用全局 source property。选择对象附加 `-fprofile-use`，其他对象普通编译，并加入
`-Werror=missing-profile -Werror=coverage-mismatch`。generate/use 使用同一全新构建
目录、同一 GCC 14.3.1、同一宏和基础 flags；不得复用 A17 已部分失败的构建树。

### A18c：验证阶梯

1. 构建日志零画像警告，严格 binary SHA256不变。
2. 高 SNR type-3，同负载 control/candidate 两轮；逐字段零差异且收益至少2%。
3. NLO有效五点和 A/B危险点。
4. C组十点与 Yukawa 2/3/4。
5. multistepmode 0/1/2/auto。
6. 42行和 safe runner 端到端。

停止条件：无法建立一一映射、任一画像警告、数值字段变化、两轮收益低于2%，或
profile-use binary只覆盖部分selected对象。

## A19：剩余编译、链接与布局候选

只评估 A1–A17 未测试的选项。候选必须先证明本机工具可用，并解释为何不会重复已失败
的全 LTO、旧画像复用或语义插入实验。优先考虑只改变代码布局或链接消除、不会改变
浮点表达式的方式。每项单独构建，不把多个 flags 捆绑测试。

固定顺序：先测 `lld + --icf=safe`，再测干净 Release `-fno-PIE/-no-pie`；只有前两项
有证据接近2%时，才考虑 action 单文件 `-fno-math-errno` 或单独代码对齐。gold ICF
仅作为 lld 的后备。禁止 `fast-math`、关闭异常/展开表和 section shuffle。

最低实验：静态资格→高 SNR两轮配对→NLO五点。数值变化立即淘汰；工具链需要修改
主项目或系统安装、预期收益低于2%、构建日志无法证明实际启用，也立即停止。

## A20：guarded 端到端流程优化

目标不是进一步放宽近似算法，而是减少必然回退点的重复工作：

- 用已有 exact/approx 双向证据审计还能否增加严格的前置 exact 区域；不得按单个
  浮点阈值向未知邻域外推。
- 检查 wrapper 是否重复解析输入、重复启动进程或重复读取参考数据。
- 分别统计首遍接受点、首遍后回退点和前置 exact 点的墙钟成本，不能只报告近似首遍。
- 任何规则只要扩大接受面、可能漏掉 B-edge/NLO边缘，或依赖事后才知道的输出量，
  都不实施。

固定优先级：首先实现批量 runner——一个 approximate CalcGW 处理多行，但每行仍独立
执行现有 guard，风险行单独 exact-fast，最终按输入顺序合并。其次才考虑 prefilter
一次加载/分类、仅凭既有双向邻域证据扩充 exact-direct 注册表，以及合并两个 Python
helper 的启动。任一方案不得改变 guard 阈值、anchor 半径或接受面。

验收以完整端到端 workload 为准：结果必须与现有 safe runner 完全相同，且总墙钟
至少改善2%。只把成本从一个进程转移到另一个进程不算收益。

## 记录与提交

每个子轮次都在 `RESEARCH_JOURNAL_ZH.md` 记录假设、唯一变量、样本、并发、状态/
history/SNR、耗时、反例、结论、产物、commit和推送状态。失败实验同样提交，避免
日后重复。只有通过全阶梯的候选才能更新默认 wrapper。
