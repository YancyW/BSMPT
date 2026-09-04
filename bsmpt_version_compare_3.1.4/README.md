# BSMPT v3.1.4 与当前主项目对照

本目录只保存版本对照的输入、输出、统计脚本和中文报告。官方 v3.1.4 源码固定在
`upstream-v3.1.4/`，tag `v3.1.4`，commit
`0b95ac4dad5bf6eee017f0c7094140b9903054d0`；源码和构建目录不纳入父仓库。

比较基准不是目录名所暗示的纯 v3.1.8：当前主项目 `CMakeLists.txt` 声明 3.3.1，
核心上游基线为 `04cb17d1233522f3c423cbd957a8922be037241e`，本次运行使用未修改的
`install/bin/CalcGW`。父仓库 HEAD 的后续提交未再改动主项目核心 C++。

## 固定计划

1. 固定官方 v3.1.4 源码、当前 HEAD 和实际版本号。
2. 由 Luna 分别审计理论/物理变化与代码/API/数值算法变化。
3. 用相同系统 GCC 14.3.1 和当前 Conan 依赖隔离构建 v3.1.4。
4. 显式关闭 Eigen/compiler vectorization；普通构建仅作为 double-free 反例保留。
5. 对高 SNR Yukawa、broad A/B/C、NLO 内外与已知危险锚点运行两个版本。
6. 比较状态、transition history、相变温度、alpha、beta/H、声速、频率与 SNR。
7. 区分绝对误差、相对误差、弱信号放大、状态翻转和运行时间。
8. 将理论、代码、实测结果、局限和复现方法写入 `VERSION_COMPARISON_ZH.md`。

当前完整结论见 `VERSION_COMPARISON_ZH.md`。原8点数值见
`version_result_summary.tsv`，新增13点见 `expanded_result_summary.tsv`；合并比较现为
21点，状态/history 翻转仍为 0/21。`summarize_expanded.py` 按表头解析宽表，避免
固定列号误读。

## 构建注意

v3.1.4 的 `conanfile.py` 明确说明 BSMPT 与 libcmaes 的 vectorization 配置不一致会
触发 double free。仅传 `-DBSMPTUseVectorization=OFF` 时，当前 Conan toolchain 仍
注入 `-ftree-vectorize`，首个测试点复现 `double free or corruption`。正式对照树
`build-v3.1.4-novec/` 因此额外使用 `-fno-tree-vectorize -DEIGEN_DONT_VECTORIZE`。
