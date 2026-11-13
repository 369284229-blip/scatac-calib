# Table 3 → 可执行映射（scatac-calib v1.1.0）

本文把论文 **Table 3（Minimum reporting checklist，15 行）** 逐项映射到
v1.1.0 的命令、输入、输出字段，并标明哪些项由 CLI 自动产出、哪些需要外部脚本。
同时给出 **Fig 1D / Fig 6 / README 中“items 1–8”** 与 Table 3 实际 15 行的编号映射。

> 只读文档：本文件不改论文正文/图/表，只说明“论文的一项要求，用什么命令能算出来”。
> 全量复现命令见 [`reproducibility.md`](reproducibility.md)；可先跑
> [`../examples/real_fixture/run_real_fixture.sh`](../examples/real_fixture/run_real_fixture.sh)
> 验证安装在真实数据子集上可用。

## 0. 计数与编号

- Table 3 实际 **15 行**（`写作/tables/table3.tsv`）。
- 决策树（`docs/decision_tree.md`，Fig 6A 的文字版）的“最小报告顺序”是 **8 项**。
- Fig 1D、Fig 6 图注、Methods §6.10 和 README 里的 “items 1–8” 指的是**决策树 8 项**，
  不是 Table 3 的 15 行。两者不是同一个编号体系；下表给出一一映射。

**一句话**：8 项是决策树的叙事顺序，15 行是 Table 3 的完整清单；8 项映射到其中 12 行，
另 3 行（checkpoint-selection、reporting template、multiple-comparison）是 Table 3
为可复现性追加的。

### 0.1 决策树 8 项 → Table 3 15 行

| 决策树 item | 内容 | 对应 Table 3 行 |
|---|---|---|
| 1 | 数据规格（细胞数、峰数、检出率、每非零元读段数） | **10** Peak-set provenance（部分；细胞/峰总数另见 Table 1） |
| 2 | 测量上限（MAPQ / 唯一比对 / M1） | **7** M1；**8** MAPQ / unique alignment |
| 3 | 主指标（逐峰 AUC、逐细胞 AUC、AP；3 seeds mean ± s.d.） | **2** Δ per-peak；**3** Δ per-cell；**12** Seeds and dispersion；**14** Average precision |
| 4 | 零模型（同一 split、同一口径） | **1** Same-split frequency null |
| 5 | ΔAUROC + 配对 bootstrap CI + p | **2**、**3**（Δ 点估计）；**13** Paired-bootstrap CI and p |
| 6 | 群级 oracle（逐细胞轴） | **4** Group-level oracle |
| 7 | 深度匹配前后 | **6** Depth matching / stratification |
| 8 | 随机峰对照（Q3 为是时） | **5** Random-peak control |

未出现在决策树 8 项中的 3 行：**9** Checkpoint-selection rule、
**11** Reporting template、**15** Multiple-comparison control。

> 投稿前的编号建议（不在本单执行）：把 Fig 1D / Fig 6 / README 的 “items 1–8”
> 改成“decision-tree items 1–8”，或把 Table 3 明确写成“8 core + 7 extensions”。
> 只改一处而不同步其余，会继续留下 046 评审指出的编号不一致。

## 1. 逐行可执行映射（15 行）

状态图例：**auto** = v1.1 直接输出该数字；**partial** = 部分由 CLI 产出/需人工补；
**external** = v1.0 不产出，需外部脚本或静态元数据。

### 1. Same-split frequency null — **auto**

| | |
|---|---|
| 命令 | `scatac-calib zero-model --matrix m_test.npz --train-matrix m_train.npz --out out/zero`（`--preds` 可省略；频率口径下 `--train-matrix` 必填） |
| 输入 | test 标签矩阵、train 矩阵（提供 `cell_freq`；train/test 峰集可不同） |
| 输出 | `zero_model.per_peak_auc`、`zero_model.per_cell_auc`、`zero_model.cell_freq_source` |
| 判定 | `cell_freq_source` 必须是 `"train"`；若为 `"test(label-informed)"` 说明用的是 inductive 模式的退化口径或细胞数不匹配，不可与 frequency 口径混比。frequency 缺 `--train-matrix` 现在**硬错误退出**（v1.1），不再静默退化 |
| 备注 | 逐峰轴用 train 的细胞检出率打分；逐细胞轴用 test 的峰检出率打分。后者是 transductive 的（046 P0-2 / 工单 047 处理），引用时须注明 |
| inductive | `scatac-calib zero-model --matrix m_test.npz --groups cell_groups.txt --null inductive --out out/ind` 只凭标签就能算逐细胞 inductive null（5-fold cell cross-fitting）；与 `--preds` 无关 |

### 2. ΔAUROC, per-peak — **auto**

| | |
|---|---|
| 命令 | `scatac-calib delta --matrix m_test.npz --preds preds.npz --train-matrix m_train.npz --out out/delta --n-boot 1000` |
| 输出 | `delta_vs_zero.peak.delta`、`.model`、`.baseline`、`.paired.ci95`、`.paired.p` |
| 备注 | 点估计 = 模型逐峰 AUC − 零模型逐峰 AUC；CI/p 来自峰级配对 bootstrap |

### 3. ΔAUROC, per-cell — **auto**

| | |
|---|---|
| 命令 | 同第 2 行，加 `--null inductive`（一次 `delta` 同时输出两轴与 inductive 块） |
| 输出 | frequency 参考：`delta_vs_zero.cell.delta`、`.paired.ci95`、`.paired.p`；inductive 主口径：`delta_vs_inductive_null.cell.{model,baseline,delta,paired.ci95,paired.p}`，null 本身在 `inductive_null.per_cell_auc` |
| 判定 | Table 3 的 `[op-default]`：Δ < 0.03 → 模型主要在复述峰频率；这是**经验判读阈值**，不是显著性门槛 |
| 全量参考 | 4,665 × 8,661：null 0.6852978；Δcell（peak-selected 3-seed 均值）+0.0193566（复现脚本 `examples/full_reproduction/`） |

### 4. Group-level oracle — **auto**（需要 `--groups`）

| | |
|---|---|
| 命令 | `scatac-calib oracle --matrix m_test.npz --preds preds.npz --train-matrix m_train.npz --groups cell_groups.txt --out out/oracle` |
| 输出 | `oracle.per_peak_auc`、`oracle.per_cell_auc`、`oracle.source`、`oracle.n_groups`、`oracle.missing_cells` |
| 泄漏口径 | `source = "train"` 才是无泄漏上界；峰集 split（train 峰数 ≠ test 峰数）时自动退化为 `"test(label-informed)"`，是上界不是可训练结果，论文中必须注明。v3 全量就是后者（14 群） |

### 5. Random-peak control — **partial**

| | |
|---|---|
| CLI 可做 | 拿到随机峰的 test 矩阵和模型预测后，`scatac-calib audit --matrix random_test.npz --preds random_preds.npz --train-matrix random_train.npz ...` 可直接给出**同峰集**零模型与 Δ（这正是 Table 3 要求的同峰集口径） |
| CLI 不做 | 随机峰集的构造、随机峰模型训练/checkpoint 选择 |
| 外部输入 | 随机峰数据集 `/root/wheat_project/02_results/wheat/dataset_random_peaks/`（构建脚本 `项目1/worker_scripts/build_dataset_random_peaks.py`） |
| 外部结果 | `p1_benchmark/random_peaks_true/seed{0,1,2}/`、`test_results_masked.json`、`zero_model.json`、`delta_vs_v3.md` |
| 输出对照 | v3 全量：Δ(v3 − random) 逐峰 +0.0587、逐细胞 +0.0376、AP +0.0450；random 自身同峰集 Δ 逐峰 +0.0039、逐细胞 −0.0940（1 seed，无 CI） |
| 备注 | 只用 1 个 seed；Table 3 的 `[op-default]` Δ > 0.01 是经验判读，不是显著性 |

### 6. Depth matching / stratification — **auto**

| | |
|---|---|
| 命令 | `scatac-calib depth --matrix m_test.npz --preds preds.npz --out out/depth --depth-bins 10`（或 `audit`） |
| 输出 | `depth_matched.per_peak_auc`、`.per_bin_auc`、`.n_bins`、`.cell_depth_median` |
| 判定 | 保留率 = (depth_matched − 0.5) / (per_peak_auc − 0.5)；`[op-default]` < 50% → 信号主要来自深度代理。v1.0 **只输出数字，不自动算保留率**，需手算 |
| 备注 | 5 箱与 10 箱不可逐位比较；报告必须写 `n_bins` |

### 7. M1 / split-half reliability — **external**

| | |
|---|---|
| v1.0 | **不实现**；`audit` 的 `reporting_checklist.md` 固定显示 `— 缺失` |
| 外部脚本 | `/root/wheat_project/项目1/代码/minloop_m1.py`（另有 `minloop_m1_annotation.py`） |
| 输入 | 同 split 的 replicate/分半矩阵（簇内 + 全细胞），不是 `m_test.npz` 单独能算 |
| 输出 | `/root/wheat_project/项目1/结果/minloop/m1_reliability.json` |
| 口径 | cluster-weighted 0.642、all-cell 0.785、AB/AD/BD 0.792/0.776/0.786；Spearman–Brown 是派生值，须标注；chance-corrected kappa 见工单 049 |
| 备注 | estimand 是**分细胞可重复性**，不是技术重复可靠性；不要写成 technical reliability |

### 8. MAPQ / unique alignment — **external**

| | |
|---|---|
| v1.0 | **不实现**；清单固定 `— 缺失` |
| 外部脚本 | `/root/wheat_project/code/rice_mapq_pilot.sh`（read-level 统计；小麦/水稻同一套） |
| 小麦 read-level | `/root/wheat_project/02_results/rice/read_level/WHEAT_AK58.mapq.tsv`（路径在 `rice/` 下，但样本是小麦 AK58；引用时写清路径） |
| 小麦逐峰 | `/root/wheat_project/02_results/wheat/dataset_v3/mapq/q0.tsv`、`q30.tsv`（chr, start, end, peak name, 读段数, 唯一/通过数；MAPQ 阈值 0 / 30） |
| 亚基因组敏感性 | `dataset_v3/subgenome_effect.json`（含 MAPQ）与 `subgenome_effect_nomapq.json`（不含）；可直接作为**亚基因组层面**的 MAPQ 敏感性 |
| 输出 | 主口径 MAPQ ≥ 30 = 47.5%；MAPQ ≤ 1 = 45.67%。**必须同时写参考基因组、aligner、MAPQ 定义**，不得混用 IRGSP-1.0 / MSU7 |
| 未做 | 全模型“MAPQ ≥ 30 过滤后重训/重评”没有现成结果；不要伪造。需要 BAM 过滤后重算，见工单 049 回执 |

### 9. Checkpoint-selection rule — **partial**

| | |
|---|---|
| CLI 可做 | 对传入的 `--preds` 计算指标；不判断该预测来自 peak-selected 还是 cell-selected |
| CLI 不做 | checkpoint 比较与选择 |
| 外部结果 | `/root/wheat_project/02_results/wheat/p1_benchmark/ablation/checkpoint_selection.tsv`（6 行：3 seeds × 2 规则；peak_auc / cell_auc_masked / AP） |
| 要求 | Table 2 每个主行必须标 `checkpoint_selection = peak` 或 `cell`；不得在同一张表里混用未标注的规则 |

### 10. Peak-set provenance (marker vs random) — **partial / static**

| | |
|---|---|
| CLI 可做 | 报告 `n_peaks`、`n_cells`（矩阵形状） |
| CLI 不做 | 峰集来源、亚基因组计数、正样本率、与随机峰的重叠 |
| 静态元数据 | `/root/wheat_project/02_results/wheat/dataset_v3/dataset_report.json`：31,741 峰（A 9,976 / B 10,457 / D 11,308）、test 正样本率 2.2123%、`peak_rule`、`split_rule`；`peaks.bed` |
| 手工计算 | 正样本率也可直接 `labels.mean()`；每非零元读段数需要原始 counts（`ad.h5ad` 已二值化，需回原始矩阵） |
| 备注 | 这一行不是 CLI 字段，需在 Methods/Table 1 手工填写 |

### 11. Reporting template — **partial**

| | |
|---|---|
| CLI 产出 | `scatac-calib audit` 同时写 `calib_report.json`（机器可读）和 `reporting_checklist.md`（人读） |
| 限制 | v1.0 的清单是 **8 行**（zero / delta / oracle / depth / AP / M1 / MAPQ / random），**不是 Table 3 的 15 行**；M1、MAPQ、random 永远显示 `— 缺失`（除非同一目录另有结果，CLI 不自动合并） |
| 结论 | “machine-readable JSON + human checklist” 成立；但“covering items 1–8 of the decision tree” 只覆盖决策树 8 项中的前 5 项自动项。要完整覆盖 Table 3，需把外部脚本输出手工并入 JSON（或等后续 `validate` 子命令） |

### 12. Seeds and dispersion — **partial**

| | |
|---|---|
| CLI 可做 | 每个 seed 单独跑一次，输出该 seed 的点估计与 bootstrap CI |
| CLI 不做 | 3 seed 的 mean ± s.d. 聚合；JSON 里没有 `mean_sd` 字段 |
| 聚合 | 按 `reproducibility.md §3` 的代码片段，`ddof=1`；**± 是 seed 间 s.d.，不是 CI** |
| v1.0 证据 | 3 seeds：逐峰 0.76459 ± 0.00108、逐细胞 0.70465 ± 0.00170 |

### 13. Paired-bootstrap CI and p — **auto**

| | |
|---|---|
| 命令 | `delta` / `audit`，`--n-boot 1000`（默认） |
| 输出 | `delta_vs_zero.{peak,cell}.paired.ci95`（2.5/97.5 分位）、`.paired.p`、`.paired.n_boot` |
| 口径 | p 在没有任何 replicate 跨 0 时截断为 `0.0`，写作 `p < 0.001`；这是**逐比较、未校正** p |
| 限制 | 当前是 i.i.d. 峰/细胞 bootstrap；block/cluster bootstrap 见工单 048，不在 v1.0 |

### 14. Average precision — **auto**

| | |
|---|---|
| 命令 | `audit` / 任意子命令加载后计算；输出在 `ap` |
| 输出 | `ap`（全矩阵 flatten 的 average precision） |
| 备注 | 必须同时报正样本率（v3 test 2.2123%）；AP 受 prevalence 影响，不能单独横比。v3 AP = 0.0879 |

### 15. Multiple-comparison control — **external / not applied**

| | |
|---|---|
| v1.0 | 不实现 BH/FDR；所有 p 都是 per-comparison |
| 论文口径 | Table 3 把它列为 recommended practice，**本文未执行**；不得写成“使用 BH q = 0.05” |
| 如需补做 | 包外对同一张表的 Δ p 做 Benjamini–Hochberg，并注明 family；脚本不属于 scatac-calib v1.0 |

## 2. 一次性跑 Table 2/3 核心项的完整命令

```bash
# 每个 seed 一次；以 seed0 为例
scatac-calib audit \
  --matrix       /root/wheat_project/02_results/wheat/dataset_v3/m_test.npz \
  --preds        /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed0/preds.npz \
  --train-matrix /root/wheat_project/02_results/wheat/dataset_v3/m_train.npz \
  --groups       /root/wheat_project/02_results/wheat/dataset_v3/cell_groups.txt \
  --out          calib_formal/seed0 \
  --n-boot 1000 --depth-bins 10 --seed 20261112
```

加 `--null inductive` 即同一命令输出 inductive null 与 `delta_vs_inductive_null`；
主口径（v1.1 起）为 inductive，`delta_vs_zero` 保留为 transductive 对照。

已保存的参考输出：`p1_benchmark/calib_formal/seed{0,1,2}/calib_report.json`。
seed0 关键值：逐峰 0.764036、逐细胞 0.703176、AP 0.087008、零模型
0.695811 / 0.687814、Δ +0.068225 / +0.015362、oracle 0.712559 / 0.777132
（`test(label-informed)`，14 群）、深度 10 箱 0.687723。

## 3. 覆盖度小结

| 状态 | 行号 | 数量 |
|---|---|---|
| auto（CLI 直接产出） | 1, 2, 3, 4, 6, 13, 14 | 7 |
| partial（CLI + 人工/外部补） | 5, 9, 10, 11, 12 | 5 |
| external（v1.0 不产出） | 7, 8, 15 | 3 |

- **7 行全自动**，**8 行需要 CLI 之外的动作**。因此“one command reproduces every
  calibration quantity”这类表述应改为“one call for the core calibration
  quantities (null / oracle / Δ / depth); M1, MAPQ, random-peak training and
  multiple-comparison control are external”。
- Table 3 的 15 行里，M1、MAPQ、random-peak 三项是 **046 评审 R3-M4** 指出的
  “永远显示缺失”的项；本文件把它们显式标为 external，并给出脚本/JSON 路径。
- `reporting_checklist.md` 的 8 行与 Table 3 的 15 行**不是同一套编号**；
  引用时写“decision-tree items 1–8”，不要写“Table 3 items 1–8”。
