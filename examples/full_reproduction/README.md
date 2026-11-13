# Full-scale inductive-null reproduction（不进包的数据）

`scatac-calib audit --null inductive` 的方法在 `examples/real_fixture/` 上可端到端
回归，但**全量论文数字不能在 100×80 子集上复现**：inductive null 的逐细胞 AUC
对峰数与细胞数都敏感（本仓库 100 峰子集的 null 是 0.6676，不是全量的 0.6853）。

全量 test split（4,665 峰 × 8,661 细胞）的标签矩阵约 6.9 MB，三个 seed 的
预测各约 142 MB；把 430 MB 放进 wheel/sdist 不合适，因此**全量 fixture 不入包**，
用本目录脚本在项目数据上只读复现：

```bash
python examples/full_reproduction/reproduce_inductive.py \
  --data-dir /root/wheat_project/02_results/wheat/dataset_v3 \
  --preds \
    /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed0/preds.npz \
    /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed1/preds.npz \
    /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed2/preds.npz
```

期望（容差 1e-3，与工单 047 回执一致）：

| 量 | 值 |
|---|---|
| inductive null（5-fold cell cross-fitting，fold_seed 20261112） | **0.6852978** |
| Δcell，per seed（peak-selected seed 0/1/2） | +0.0178785 / +0.0212185 / +0.0189728 |
| Δcell，3-seed 均值 | **+0.0193566** |

脚本退出码 0 表示两项都在容差内。数据只读，不写项目目录（`--out` 可选）。

> 单次 CLI 运行只对应一个 `--preds`，因此输出的是单个 seed 的 Δ；+0.0194 是
> 3 个 seed 的均值，按 `docs/reproducibility.md` 的“每个 seed 跑一次、包外聚合”
> 口径得到。null 与 seed 无关。
