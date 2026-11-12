# 从原始矩阵到 Table 2 的最短路径

本文给出用 `scatac-calib` v1.0.0 复现论文 **Table 2**（主指标 / 零模型 /
oracle / Δ 与 CI）与 **Table 3** 清单项的最小步骤。所有命令只读取已有
矩阵和预测，不训练、不重抽峰集。

## 0. 前置物

每个分析的 split 需要：

| 文件 | 形状 / 内容 | 来源 |
|---|---|---|
| `m_test.npz` | 峰 × 细胞，CSR，二值或 counts | 数据准备阶段已固定，测试集**不重抽** |
| `m_train.npz` | 峰 × 细胞，CSR | 同一 split 的训练集；峰集可与 test 不同 |
| `cell_groups.txt` | 每行 `barcode<TAB>cluster`，行数 = test 细胞数 | Leiden/细胞类型标注 |
| `preds.npz`（每 seed 一个） | 峰 × 细胞，键 `preds` | 已训练模型的 test 预测；不重跑训练 |

> 关键约束：`m_test.npz`、`cell_groups.txt` 的行列顺序必须严格对应。
> `preds` 与 `m_test` 形状必须一致，CLI 会检查并报 `shape mismatch`。

## 1. 安装（一次）

```bash
python -m venv --system-site-packages /tmp/calib_venv
/tmp/calib_venv/bin/pip install --no-build-isolation --no-deps -e /path/to/scatac_calib
/tmp/calib_venv/bin/scatac-calib --help
```

依赖 `numpy` / `scipy` / `scikit-learn` 已在环境里；若没有，去掉
`--system-site-packages` 与 `--no-deps` 正常安装即可。

## 2. 每个 seed 跑一次 audit

```bash
for s in 0 1 2; do
  scatac-calib audit \
    --matrix       /path/to/dataset_v3/m_test.npz \
    --preds        /path/to/p1_benchmark/bootstrap/seed${s}/preds.npz \
    --train-matrix /path/to/dataset_v3/m_train.npz \
    --groups       /path/to/dataset_v3/cell_groups.txt \
    --out          /path/to/calib_formal/seed${s} \
    --n-boot 1000 --depth-bins 10 --seed 20261112
done
```

每个 `seed${s}/` 产出：

- `calib_report.json` —— 下表所有字段
- `reporting_checklist.md` —— Table 3 的 8 项勾选

参考耗时：4,665 峰 × 8,661 细胞、1,000 bootstrap、10 深度箱约 7–8 分钟/seed（单核）。

## 3. 字段到 Table 2 的映射

| Table 2 单元 | `calib_report.json` 字段 |
|---|---|
| 逐峰 AUROC | `metrics.per_peak_auc` |
| 逐细胞 AUROC | `metrics.per_cell_auc` |
| AP | `ap` |
| 零模型 逐峰 / 逐细胞 | `zero_model.per_peak_auc` / `zero_model.per_cell_auc` |
| 零模型来源（train / test） | `zero_model.cell_freq_source` |
| Δ逐峰 / Δ逐细胞 + CI + p | `delta_vs_zero.peak.delta` / `.paired.ci95` / `.paired.p`；cell 同构 |
| 均值 CI（非 Δ） | `ci95.per_peak_auc` / `ci95.per_cell_auc` |
| oracle 逐峰 / 逐细胞 | `oracle.per_peak_auc` / `oracle.per_cell_auc` |
| oracle 泄漏口径 | `oracle.source`（`train` = 无泄漏；`test(label-informed)` = 上界） |
| 深度分层逐峰 | `depth_matched.per_peak_auc`、`.per_bin_auc`、`.n_bins` |

**Table 2 的“均值 ± s.d.”不来自单次运行**：工具输出每个 seed 的点估计，
3 个 seed 的均值与样本标准差（`ddof=1`）在包外聚合：

```python
import json
from pathlib import Path
import numpy as np

rows = [json.loads(Path(f"calib_formal/seed{s}/calib_report.json").read_text())
        for s in (0, 1, 2)]
for axis in ("peak", "cell"):
    vals = [r["metrics"][f"per_{axis}_auc"] for r in rows]
    print(axis, f"{np.mean(vals):.5f} ± {np.std(vals, ddof=1):.5f}")
```

注意：**CI 只用 seed-0 预测**（`--n-boot 1000`），不覆盖训练随机性；
3 个 seed 不合成合并 CI。Table 2 的 ± 是 seed 间 s.d.，不是 CI。

## 4. 口径检查清单（复跑前逐条确认）

1. **checkpoint 选择**：peak-selected 与 cell-selected 是两套预测，不能混用。
   Table 2 主行是 cell-selected；**全部 Δ/CI 只在 peak-selected 上报告**。
   传入的 `preds.npz` 必须与你要报告的行一致。
2. **同 split**：零模型、oracle、Δ 必须和主指标用同一个 `m_test.npz` 与同一个
   `cell_groups.txt`；不要为不同行重抽测试集。
3. **零模型来源**：峰集 split 时 `cell_freq_source` 必须是 `"train"`。
   若显示 `"test"`，说明没传 `--train-matrix` 或细胞数不匹配，数字不可比。
4. **oracle 泄漏**：`oracle.source == "train"` 才是可训练意义上的上界；
   若为 `"test(label-informed)"`，论文中必须注明。
5. **深度箱数**：5 箱与 10 箱不可逐位比较；报告时写明 `n_bins`。
6. **bootstrap RNG**：固定 `--seed`。与另一脚本对比时，点估计应逐位一致，
   CI 可能差 1e-4 量级，这是 RNG 流不同，不是口径冲突。

## 5. 回归自检（v1.0.0 已验证）

031 正式跑（3 seeds，1000 bootstrap，10 深度箱）与
`bootstrap_summary.json` 的对比判据：

- Δ 点估计差 ≤ 1e-4 → 实测 **3 seeds × 2 轴全部为 0**
- CI 端点差 ≤ 1e-3 → 实测最大 **1.51e-4**
- 逐峰/逐细胞 AUC 与 `B_summary.md` 差 < 1e-7
- 零模型 / oracle 差 0

期望值（与 `回执_031.md` 一致）：

| 量 | 值 |
|---|---|
| 逐峰 AUC（3 seed 均值 ± s.d.） | 0.76459 ± 0.00108 |
| 逐细胞 AUC | 0.70465 ± 0.00170 |
| Δ逐峰 | +0.06878 ± 0.00108 |
| Δ逐细胞 | +0.01684 ± 0.00170 |
| 零模型 逐峰 / 逐细胞 | 0.69581 / 0.68781 |
| oracle 逐峰 / 逐细胞 | 0.71256 / 0.77713 |
| 深度 10 箱保留率 | 71.10% / 71.75% / 71.01% |

任一 seed 复跑后若点估计偏离上表，先查 checkpoint 选择与 split，再查输入
行列顺序；不要通过改算法或改 `--seed` 去凑数字。

## 6. Table 3 中工具不能自动产出的项

清单里以下项在 v1.0 显示“— 缺失”，需外部脚本：

- **M1 分半信度**：需要同 split 的 replicate 矩阵（簇内 + 全细胞，Spearman–Brown 校正）。
- **MAPQ≤1 / 唯一比对比例**：需要 BAM 与 read-level 统计。
- **随机峰对照**：需要按同 split 规则重抽峰集并重新得到预测；本工具只评测给定矩阵。

这些项不阻塞 Table 2 的复现，但 Table 3 的完整勾选需要它们。
