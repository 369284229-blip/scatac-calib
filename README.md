# scatac-calib v1.0.0

scATAC-seq **序列模型校准评测** 的命令行工具。回答一个问题：

> 模型报出的 AUROC，有多少是数据本身的白送（峰频率、细胞深度），有多少是序列信号？

一条命令给出：零模型、群级 oracle、ΔAUROC + 配对 bootstrap CI、深度分层 AUC，
以及一份最低报告清单。算法行为自 v0.1 起未改变；v1.0 只是可安装包 + 测试 + demo + 文档。

- Python ≥ 3.9，依赖 `numpy` / `scipy` / `scikit-learn`
- 五个子命令：`audit` / `zero-model` / `oracle` / `delta` / `depth`
- 许可证 MIT，引用信息见 `CITATION.cff`

## 安装

从源码目录（推荐开发时用 editable）：

```bash
cd scatac_calib
python -m pip install -e .
scatac-calib --help
```

从 wheel：

```bash
python -m pip install dist/scatac_calib-1.0.0-py3-none-any.whl
scatac-calib --help
```

不安装也可以直接跑（需在源码根目录）：

```bash
python -m scatac_calib --help
```

可选：读 `.h5ad` 需要 `scanpy`（`pip install "scatac-calib[h5ad]"`）。只读 `.npz`
时不需要。

## 输入格式

| 参数 | 说明 |
|---|---|
| `--matrix` | test 标签矩阵，**峰 × 细胞**；`.npz`（CSR 四件套 `data`/`indices`/`indptr`/`shape`）或 `.h5ad` |
| `--preds` | 模型预测，与 `--matrix` 同形状；`.npz`（默认键 `preds`）或 `.npy` |
| `--preds-key` | `.npz` 中的预测键名，默认 `preds` |
| `--train-matrix` | 训练矩阵（峰 × 细胞），用于零模型的细胞检出率；峰集与 test 不同也可以 |
| `--groups` | test 细胞分组，`barcode<TAB>cluster`，每行一条；用于群级 oracle |
| `--train-groups` | 训练细胞分组；缺省时复用 `--groups` |
| `--out` | 输出目录 |
| `--n-boot` | bootstrap 次数，默认 1000 |
| `--depth-bins` | 深度分箱数，默认 10 |
| `--seed` | bootstrap RNG 种子，默认 20261112 |

写成 CSR `.npz`：

```python
import numpy as np
from scipy import sparse
m = sparse.csr_matrix(Y)          # Y: 峰 x 细胞
np.savez_compressed("m_test.npz", data=m.data, indices=m.indices,
                    indptr=m.indptr, shape=m.shape)
np.savez_compressed("preds.npz", preds=P)
```

> 峰集 split（train/test 峰不同）时，零模型的 `cell_freq` 只能来自训练矩阵列均值；
> `peak_freq` 来自 test 峰检出率。这与项目 S-013/S-014 的口径完全一致。

## CLI 用法

一次跑完全部（推荐）：

```bash
scatac-calib audit \
  --matrix m_test.npz \
  --preds preds.npz \
  --train-matrix m_train.npz \
  --groups cell_groups.txt \
  --out out/wheat_seed0 \
  --n-boot 1000 --depth-bins 10
```

输出：

- `calib_report.json`：全部数字
- `reporting_checklist.md`：最低报告清单（✅ / — 缺失）

单跑某一项（参数相同）：

```bash
scatac-calib zero-model --matrix m_test.npz --preds preds.npz \
  --train-matrix m_train.npz --out out/zero
scatac-calib oracle    --matrix m_test.npz --preds preds.npz \
  --train-matrix m_train.npz --groups cell_groups.txt --out out/oracle
scatac-calib delta     --matrix m_test.npz --preds preds.npz \
  --train-matrix m_train.npz --out out/delta
scatac-calib depth     --matrix m_test.npz --preds preds.npz \
  --out out/depth --depth-bins 10
```

## Demo（一条命令）

包内有两个 100 峰 × 80 细胞输入，用途不同，**不要混用**：

| 目录 | 数据来源 | 训练矩阵 | 分组 | 用途 |
|---|---|---|---|---|
| `examples/demo/` | **synthetic demo**：固定种子合成（4 群，25 个群 marker 峰/群） | 合成，同一生成过程 | 4 个合成群 | 零依赖教程，验证流程 |
| `examples/real_fixture/` | **real-data fixture**：小麦 `dataset_v3` test 的固定索引子集（100 峰，chr3） | 真实，23,268 个 train 峰（与 test 峰按染色体不相交） | 真实 14 个 Leiden 群 | 在真实数据上核对安装与输出 |

### synthetic demo

不依赖任何真实数据即可复现完整流程：

```bash
cd examples/demo
./run_demo.sh            # 默认输出 demo_out/
```

等价命令：

```bash
scatac-calib audit \
  --matrix m_demo.npz --preds preds_demo.npz \
  --train-matrix m_train_demo.npz --groups groups_demo.txt \
  --out demo_out --n-boot 200 --depth-bins 5
```

参考结果与容差见 [`examples/demo/expected.md`](examples/demo/expected.md)。
一句话摘要：逐峰 0.62766、逐细胞 0.74127；零模型 0.50583 / 0.68782；
Δ +0.12183 / +0.05345（p < 0.001）；oracle 逐细胞 0.75501；深度 5 箱逐峰 0.66879。

数据可确定性地重建：

```bash
python examples/demo/make_demo_data.py --seed 20261112
```

### real-data fixture

`examples/real_fixture/` 是真实数据的固定子集：从 `dataset_v3` test（4,665 峰）
按 `np.linspace` 取 100 峰、从 14 个 Leiden 群各取 5–6 个细胞（共 80），
预测取 seed-0 peak-selected `preds.npz` 的对应子块。行索引、选择规则与
SHA256 见 `selection_manifest.json`。

```bash
cd examples/real_fixture
./run_real_fixture.sh    # 默认输出 fixture_out/
```

参考数字、容差与“为什么这个小子集不能当科学结论”见
[`examples/real_fixture/expected_real.md`](examples/real_fixture/expected_real.md)。
该 fixture 的训练峰与 test 峰按染色体不相交，因此 oracle 显示
`test(label-informed)`（上界，不是可训练结果）；这与全量 v3 的口径一致。

## 测试（不依赖 pytest）

```bash
cd scatac_calib
pip install -e .
python -m unittest discover -s tests -v
# 或
./tests/run_tests.sh
```

- `tests/test_zero_model.py`：200 峰 × 400 细胞随机二值矩阵上，零模型逐峰/逐细胞
  AUC 落在 0.45–0.55；模型等于零模型时 Δ = 0；全零/全一峰与细胞被跳过（NaN）。
  用 400 而非 100 细胞是因为 `peak_freq` 来自 test 矩阵本身：100 细胞时逐细胞轴的
  有限样本自评分偏差约 +0.056，400 细胞时 < 0.03，0.45–0.55 才是有效判据。
- `tests/test_cli_smoke.py`：在 demo 数据上端到端跑 `audit`，检查 JSON/清单字段。
- `tests/test_real_fixture.py`：在 `examples/real_fixture/` 的真实数据子集上端到端
  跑 `audit`，点估计与 `expected_real.md` 的容差为 1e-6（`--n-boot 50` 只影响 CI，
  不影响点估计）。共 10 个测试。

干净 venv 验证命令：

```bash
python -m venv --system-site-packages /tmp/calib_venv
/tmp/calib_venv/bin/pip install --no-build-isolation --no-deps -e .
/tmp/calib_venv/bin/python -m unittest discover -s tests -v
```

## 与论文 Table 2 / Table 3 的对应

| 论文表 | 内容 | 工具字段 |
|---|---|---|
| Table 2（主指标） | 逐峰 / 逐细胞 AUROC、AP（3 seeds 均值 ± s.d.） | `metrics.per_peak_auc`、`metrics.per_cell_auc`、`ap`；每个 seed 跑一次后在包外聚合 |
| Table 2（零模型行） | 同 split、同峰集的频率零模型 | `zero_model.per_peak_auc`、`zero_model.per_cell_auc`、`zero_model.cell_freq_source` |
| Table 2（oracle 行） | 群级 oracle（逐峰/逐细胞） | `oracle.per_peak_auc`、`oracle.per_cell_auc`、`oracle.source`、`oracle.n_groups` |
| Table 2（Δ 与 CI） | ΔAUROC + 配对 bootstrap 95% CI + p | `delta_vs_zero.{peak,cell}.delta`、`.paired.ci95`、`.paired.p` |
| Table 3（最低报告清单） | 8 项报告要求 | `reporting_checklist.md`；`depth_matched.*` 填深度项 |
| Table 3（M1 / MAPQ / random peaks） | 需要 replicate 矩阵 / BAM / 重抽峰集 | v1.0 未实现，清单中显示“— 缺失”，需外部脚本补齐 |

最短复现路径（原始矩阵 → Table 2）见 [`docs/reproducibility.md`](docs/reproducibility.md)。

## 结果口径

| 项 | 定义 |
|---|---|
| 逐峰 AUC | 对每个峰，用模型分数给细胞排序，算跨细胞 AUROC；对峰取平均 |
| 逐细胞 AUC | 对每个细胞，用模型分数给峰排序，算跨峰 AUROC；对细胞取平均 |
| 零模型 | 逐峰轴用**细胞检出率**打分；逐细胞轴用**峰检出率**打分 |
| 群级 oracle | 用每个细胞所属群的峰均值打分（`source=train` 无泄漏；峰集不同时退化 `test(label-informed)`，须在文中注明） |
| 深度分层 AUC | 按细胞深度分位分箱，箱内算逐峰 AUC 再平均 |
| 配对检验 | `--n-boot` 次配对 bootstrap，双侧 p |

## 已验证（031 回归，2026-11-12）

在 L-wdl `bioinfo_env` 上对小麦数据跑正式版（1000 bootstrap / 10 深度箱），3 个 seed：

| seed | 逐峰 AUC | 逐细胞 AUC | Δ逐峰 | Δ逐细胞 |
|---|---|---|---|---|
| 0 | 0.76404 | 0.70318 | +0.06823 | +0.01536 |
| 1 | 0.76584 | 0.70652 | +0.07002 | +0.01870 |
| 2 | 0.76389 | 0.70427 | +0.06808 | +0.01646 |
| 均值 ± s.d. | 0.76459 ± 0.00108 | 0.70465 ± 0.00170 | +0.06878 ± 0.00108 | +0.01684 ± 0.00170 |

- 所有 Δ 的配对 bootstrap `p = 0.0`（p < 0.001）。
- 零模型（同 split / 同峰集）：逐峰 **0.69581**、逐细胞 **0.68781**（三 seed 相同）。
- 群级 oracle（14 群，test label-informed）：逐峰 **0.71256**、逐细胞 **0.77713**。
- 深度 10 箱保留率（相对 0.5 的增益）：71.10% / 71.75% / 71.01%（均值 71.29%）。
- **回归判据**：Δ 点估计差 ≤ 1e-4、CI 端点差 ≤ 0.001。实测 3 seeds × 2 轴的 Δ 点估计差
  **全部为 0**；CI 端点最大差 **1.51e-4**（seed0 Δpeak 上端）→ **通过**，无算法漂移。
- 与 `bootstrap_summary.json` 的差异仅来自 bootstrap RNG 流（不同脚本），点估计逐位一致。
- README 旧版 5 箱数字（0.69441 / 73.4%）与 031 的 10 箱（0.68772 / 71.3%）分箱数不同，
  同量级，不逐位比较。

## 边界

- v1.0 不含 M1 分半信度、MAPQ 统计与随机峰对照（需要 replicate 矩阵 / BAM / 重抽峰集，
  接口已在清单中预留）。
- 一次只读一个 test 矩阵；超大矩阵请先降采样或分块。
- 输出数字只描述**给定 split 上**的校准情况，不构成对模型架构的评价。
