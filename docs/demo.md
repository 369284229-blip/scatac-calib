# 30 分钟 demo

目标：别人拿到自己的 cell × peak 矩阵和模型预测，30 分钟内跑出零模型、ΔAUROC、深度控制。

## 0–5 min 准备输入

三样东西：

1. **test 矩阵**：峰 × 细胞，`.npz`（CSR）或 `.h5ad`
2. **预测矩阵**：与 test 同形状，`.npz` 里键名 `preds`
3. **训练矩阵**：用于零模型的细胞检出率（峰集与 test 不同也没关系）

可选：细胞分组文件（群级 oracle 用）。

CSR 写法参考：

```python
import numpy as np
from scipy import sparse
m = sparse.csr_matrix(Y)          # Y: 峰 x 细胞，0/1 或 counts
np.savez_compressed("m_test.npz", data=m.data, indices=m.indices,
                    indptr=m.indptr, shape=m.shape)
```

## 5–10 min 安装与自检

```bash
cd scatac_calib
python -m scatac_calib --help
```

## 10–25 min 跑

```bash
python -m scatac_calib audit \
  --matrix m_test.npz \
  --preds preds.npz \
  --train-matrix m_train.npz \
  --groups my_groups.txt \
  --out out/my_run \
  --n-boot 1000 --depth-bins 10
```

参考耗时：4,665 峰 × 8,661 细胞、10 深度箱约 **7–8 分钟**（单核）。

## 25–30 min 读结果

打开 `out/my_run/calib_report.json`：

| 字段 | 怎么读 |
|---|---|
| `metrics.per_peak_auc` | 主指标（逐峰轴） |
| `zero_model.per_peak_auc` | 白送部分；两者差 = ΔAUROC |
| `delta_vs_zero.peak.paired` | 95% CI 与 p；CI 不含 0 才算显著 |
| `oracle.per_cell_auc` | 已知细胞群的上限 |
| `depth_matched.per_peak_auc` | 去掉深度效应后还剩多少 |

再看 `reporting_checklist.md`，把缺失项补齐（M1、MAPQ、随机峰对照）。

## 常见坑

- **形状反了**：矩阵必须是 峰 × 细胞；CLI 会报 shape mismatch。
- **预测没有 sigmoid**：本工具只做排序，logits 也能算 AUC，但不要混入 NaN。
- **train/test 峰集不同**：正常。零模型会自动用 train 的细胞检出率。
- **群级 oracle 显示 `test(label-informed)`**：说明峰集不重合，oracle 用 test 群均值算，
  是上界不是可训练结果，论文里必须注明。
