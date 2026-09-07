# Region-v4 认证域优化提交报告

日期：2026-09-07

## 提交目标

在不修改主项目和冻结严格路径的前提下，提交一套可执行的非严格路由方案。首要契约
是严格有限`SNR>0`当且仅当候选有限`SNR>0`；每个认证域独立要求FP=FN=0。

## 当前方案

- 已知bounce/数值危险域输入级直接严格。
- 所有点先运行严格NLO-only前缀，消除近似NLO边界FP。
- 正点必须有三温shadow证书并命中显式正点认证域。
- “中心失败、两邻温成功”只在B3/B4两个认证点开放。
- 近似`no_coex_pair`只有命中显式coexistence认证域才接受。
- 域外、诊断缺失/分歧、晚期失败全部严格回退。

## 多样性扩展

本轮新增四个彼此分离的局部区域，每区6点：

| 区域 | 类型 | 严格/候选 | FP/FN | v4修正降时 |
|---|---|---:|---:|---:|
| E14 broad-B row1 | coexistence失败 | 6 TN | 0/0 | 14.00% |
| E15 broad-C row1 | coexistence失败 | 6 TN | 0/0 | 19.18% |
| E16 broad-A row1 | 正点 | 6 TP | 0/0 | 60.22% |
| E17 broad-C row5 | 正点 | 6 TP | 0/0 | 61.42% |

修正降时包含从诊断侧栏补回的严格NLO-only前缀。四组均先经strict/v3探索，再写入
认证表，最后使用正式v4逐点实际回归；不是仅做离线距离匹配。

## 累计证据

- 去重严格配对192点：TP=68、TN=124、FP=0、FN=0。
- 核心42点已有路由：完整严格4/42，v3同计算路径实测约降低45%。
- 192点并非独立同分布随机样本，其中A4 coexistence局部加密较多；因此只能证明
  已采样认证域，不代表全参数空间零错误率。
- classified数据经THDMTools预筛，只用于发现候选/反例，从未成为生产路由特征。

## 入口

```bash
bsmpt_speed_lab/run_calcgw_approx_region_v4_certified.sh \
  --model=r2hdm --input=input.tsv --output=output.tsv \
  --firstline=2 --lastline=2
```

外层批量运行必须使用`parallel_calcgw.py --jobs 2`或更低。runner通过输出锁和原子
替换防止迟到进程覆盖reference；比较器检查输入一致性、transition schema和文件哈希。
