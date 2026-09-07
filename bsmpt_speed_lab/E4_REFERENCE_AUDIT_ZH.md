# E4 strict/v2/v3 参考一致性审计

日期：2026-09-06。本文只检查已有文件和脚本，不运行 CalcGW；不修改主项目。

## 结论

“v2 为 0/10 个 SNR 幅值误差超过 10%，v3 为 9/10”首先是参考文件身份漂移造成的，不能解释为 v3 算法使 E4 幅值突然恶化。

同时，E4 的 strict 结果本身相对 v2/v3 出现了明显的 transition/SNR 结构变化，因此该批次也暴露出 strict 重复结果不稳定（或至少未证明可重复）。现有证据不能证明 legacy shadow 改变了状态；相反，v2 与 v3 的非 runtime 输出逐字段相同，shadow 只改变诊断和耗时。

## 可复核证据

### 1. 文件写入顺序不支持把旧 v2 摘要和当前 strict 配对

| 文件 | mtime（+0800） |
|---|---|
| 输入 `e4_broad_c_row3_delta1e4_10.tsv` | 2026-09-05 21:08:35 |
| v2 输出 | 2026-09-06 08:39:07 |
| strict 输出 | 2026-09-06 08:47:21 |
| v3 输出 | 2026-09-06 09:23:50 |
| v3 评估明细/摘要 | 2026-09-06 09:24:38 |
| v2 评估明细/摘要（后来重写） | 2026-09-06 20:11:31 |

旧 v2 结果被保存在 `region_v2_total_72_eval_details.tsv` 的 E4 行中，10 个值为
`0.0005106, 0.0012731, ..., 0.0021987`，即全部低于 10%。该聚合文件与当前
`e4_*_region_v2_eval_details.tsv` 不一致：后者已经变成了 9 个超 10% 的 v3 式值。

当前两个局部评估明细文件 SHA-256 完全相同：

```text
e3b7d0d7e0b4ac855f017ef0d3145f9c9bc738af75e0430bdd351dbca8b2dc90  e4_*_region_v2_eval_details.tsv
e3b7d0d7e0b4ac855f017ef0d3145f9c9bc738af75e0430bdd351dbca8b2dc90  e4_*_region_v3_eval_details.tsv
```

当前 v2 局部摘要也已是 `secondary_snr_amplitude_over_10pct: 9`，而历史 v2 总汇总的 E4 行仍是 0/10。这是评估产物被同名输出重写的直接证据，足以解释“v2 0、v3 9”的表面矛盾。

### 2. v2 和 v3 实际计算结果相同（除 runtime）

逐行比较 `e4_broad_c_row3_delta1e4_10_region_v2.tsv` 与 v3 文件的全部公共字段，排除 `runtime` 后 10/10 行差异数均为 0。两者每行都有
`transition_history=0-(0)->1`，SNR 总值约 `3.78–3.82e-12`。

因此 v3 三温 shadow 没有造成 E4 的结果状态或 SNR 改变；它只在诊断中记录三温成功，并带来约 0.1–0.5 s 的 runtime 差异。这个结论也与 wrapper 一致：v3 只是通过环境变量替换 runner/decision，最终输出仍由同一 v2/v1 路径产生。

### 3. 当前 strict reference 与 v2/v3 不是同一物理分支的数值形态

当前 strict 的 10 行中，除第 2 行外大多为
`transition_history=0-(0)->1-(1)->2`，第 2 行为 `0`；v2/v3 全部为
`0-(0)->1`。strict `SNR(LISA-3yrs)_0` 为约 `1.5e-37–1.7e-34`（第 2 行为 `nan`），而 v2/v3 为约 `3.8e-12`。

对同一 `_0` 字段直接相除，9 个可比较行的倍率约为 `2.3e22–2.5e25`；这不是正常的浮点舍入。评估脚本 `compare_qualitative_outcomes.py` 只把“任意 transition 有正 SNR”作为 qualitative positive，并按同名 SNR 字段计算次级误差；它不会要求 `transition_history` 或 transition 数量一致。因此 strict 仍被标成 positive，但会把两种不同 transition 结构误报成幅值差异。第 2 行 strict 为 `nan`，被脚本跳过，故恰好得到 9/10 而不是 10/10。

### 4. 是否是 shadow 对状态有副作用？

现有证据不支持该解释：

* v2 与 v3 输出的全部非 runtime 字段逐字段相同，包含所有 status、transition history 和 SNR；
* v3 的多温诊断只写到 stderr 临时日志，decision 脚本据此决定接受/回退，不修改 strict 文件；
* `run_calcgw_exact_fast.sh` 不设置 shadow 环境变量，strict 输出在 v2 之后才写出。

所以可以确认的是“同名参考/评估文件发生了迟到批次覆盖”；strict 数值差异则应标记为“不同 strict run/未验证可重复性”，不能归因给 shadow。仅凭 mtime 无法证明两个 CalcGW 进程是否曾经同时运行，但没有必要依赖并发假设：迟到 strict 写入和同名评估重写已经足以破坏配对。

## 修复方案

1. **每次运行使用不可复用的 run id 目录/前缀。** 输出至少包含 `e4.../strict/<run-id>.tsv`、`.../v2/<run-id>.tsv`、`.../v3/<run-id>.tsv`；评估明细不得再次写入固定的 `e4_*_eval_*` 名称。禁止 wrapper 和评估器覆盖已有文件；存在则失败。
2. **评估前写 manifest。** 记录输入、strict、approx、diagnostics 的绝对路径、SHA-256、mtime、命令行、binary SHA-256、环境变量和开始/结束时间；评估摘要复制这些 hash。配对时只接受 manifest 中的 hash，不按文件名或行号猜测。
3. **strict reference 做完整性门禁。** 评估前检查输入参数列逐行相等；检查 transition id/history 集合以及每个 transition 的 status/SNR schema。strict 与 approx 的 transition 集合不同就输出 `reference_structure_mismatch`，禁止产生“幅值精度”统计，必须先重跑成同一批次的 strict/approx 配对。
4. **对 strict 重复做独立复现。** 同一输入连续运行至少两次，保存各自独立输出；若 transition/history 或 SNR 超过预设容差，标记该点为 strict-unstable，不能用作 v2/v3 幅值证书。不要用当前 E4 的 strict 文件覆盖旧 reference。
5. **让评估脚本原子且不可覆盖。** `compare_qualitative_outcomes.py` 应使用 `O_EXCL`/临时文件后原子 rename，并在目标已存在时退出；摘要中增加 `reference_sha256`、`candidate_sha256`、`structure_mismatch_count` 和 `evaluation_run_id`。

在完成上述身份保护和 strict 重复验证前，E4 只能报告：v2/v3 的 qualitative 结果均为 10/10 positive；E4 的 SNR 幅值比较目前不可用。不能把当前 9/10 作为 v3 幅值失败，也不能把旧 v2 的 0/10 作为与当前 strict reference 的有效配对。
