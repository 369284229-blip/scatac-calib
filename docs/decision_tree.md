# 校准评测决策树（Fig 6A 文字版）

> 给任何做 scATAC 序列→可及性预测的人。照着走，就不会把“峰频率白送”当成模型能力。

## Q1. 你的矩阵是否稀疏（每非零元 ≈ 1 条 read，检出率 <5%）？

- **是** → 先做**测量上限**描述：每非零元读段数、唯一比对比例、split-half 一致性
  - 有多倍体/同源组 → 追加 **M1 分半信度**（簇内 + 全细胞），报告 Spearman–Brown 校正值
  - 有 BAM → 报告 **MAPQ≤1 比例**与唯一比对比例
- **否** → 跳过 M1，直接 Q2

## Q2. 主指标是逐峰还是逐细胞？

| 主指标 | 必须报告 | 原因 |
|---|---|---|
| **逐峰 AUC** | 零模型、ΔAUROC、**深度分层 AUC** | 深细胞更容易被判对；不分层会虚高 |
| **逐细胞 AUC** | 零模型（**inductive** 优先）、ΔAUROC、**群级 oracle** | 峰频率可白送大部分 AUC；oracle 给出“知道细胞类型”的上限 |
| 两者都报 | 上面全部 + AP | 审稿人一定会问另一个轴 |

## Q3. 你的峰集是否用到了注释/差异分析来选峰？

- **是** → 必须加 **随机峰对照**（同数量、同 split 规则重抽），报告 ΔAUROC(v3 − random)
  - Δ > 0.01 → 正文写明“marker 选峰抬高了指标”
- **否** → 可选

## Q4. 是否有第二个物种/条件要做对照？

- **是** → **深度匹配**后再比较：按细胞深度分位匹配，报告匹配前后两套数字
  - 深度不同还会污染“零模型 AUROC”本身（水稻 pilot 0.83 vs 小麦 0.69 就是例子）
- **否** → 略

## 最小报告顺序（Table 3）

1. 数据规格：细胞数、峰数、检出率、每非零元读段数
2. 测量上限：MapQ / 唯一比对 / M1（若适用）
3. 主指标：逐峰 AUC、逐细胞 AUC、AP（3 seeds mean±sd）
4. 零模型：同一 split、同一口径
5. ΔAUROC + 配对 bootstrap CI + p
6. 群级 oracle（逐细胞轴）
7. 深度匹配前后
8. 随机峰对照（若 Q3 为是）

## 判读规则

- **ΔAUROC(逐细胞) < 0.03**：模型基本只在复述峰频率，不能宣称“预测了细胞状态”
- **深度匹配后保留率 < 50%**（相对 0.5 的增益）：序列信号主要来自深度代理，必须写明
- **oracle(逐细胞) > 模型**：仍存在“知道细胞类型就能拿到”的信号，模型没有利用
- **随机峰对照 Δ > 0.01**：峰集选择效应显著，必须降级主结论

## 工具

```bash
# frequency null（逐细胞轴 transductive，作对照）
python -m scatac_calib audit --matrix ... --preds ... --train-matrix ... --groups ...

# inductive null（5-fold cell cross-fitting，逐细胞主口径）
python -m scatac_calib audit --matrix ... --preds ... --train-matrix ... --groups ... \
  --null inductive
```

> 逐细胞 Δ 主口径用 `--null inductive`：折内细胞的 null 分数只由其余折细胞估计，
> 避免“null 已经知道 test labels”的循环论证。frequency 口径保留作对照。

输出 `calib_report.json` + `reporting_checklist.md`，直接对应上面 1–8 项。
