# 最低报告清单：m_demo.npz

| 项目 | 说明 | 状态 |
|---|---|---|
| zero_model | 零模型 AUROC（逐峰/逐细胞） | ✅ |
| delta_vs_zero | ΔAUROC 与配对 bootstrap CI/p | ✅ |
| oracle | 群级 oracle（逐峰/逐细胞） | ✅ |
| depth_matched | 深度分层后的逐峰 AUROC | ✅ |
| ap | 平均精度 AP | ✅ |
| m1 | M1 分半信度（有同源组时） | — 缺失 |
| mapq | MAPQ≤1 比例 / 唯一比对比例（有 BAM 时） | — 缺失 |
| random_peaks | 随机峰对照（可选） | — 缺失 |

## 主要数字

- **n_peaks**: 100
- **n_cells**: 80
- **metrics.per_peak_auc**: 0.6276565763297718
- **metrics.n_valid_peaks**: 100
- **metrics.per_cell_auc**: 0.7412674975299385
- **metrics.n_valid_cells**: 80
- **ap**: 0.4471357130693913
- **zero_model.per_peak_auc**: 0.5058272394321816
- **zero_model.per_cell_auc**: 0.6878211592779424
- **delta_vs_zero.peak.model**: 0.6276565763297718
- **delta_vs_zero.peak.baseline**: 0.5058272394321816
- **delta_vs_zero.peak.delta**: 0.12182933689759023
- **delta_vs_zero.peak.paired.p**: 0.0
- **delta_vs_zero.peak.paired.mean**: 0.12119307568543164
- **delta_vs_zero.peak.paired.n**: 100
- **delta_vs_zero.peak.paired.n_boot**: 200
- **delta_vs_zero.cell.model**: 0.7412674975299385
- **delta_vs_zero.cell.baseline**: 0.6878211592779424
- **delta_vs_zero.cell.delta**: 0.05344633825199607
- **delta_vs_zero.cell.paired.p**: 0.0
- **delta_vs_zero.cell.paired.mean**: 0.053601999975413046
- **delta_vs_zero.cell.paired.n**: 80
- **delta_vs_zero.cell.paired.n_boot**: 200
- **oracle.n_groups**: 4
- **oracle.missing_cells**: 0
- **oracle.per_peak_auc**: 0.6465986682226569
- **oracle.per_cell_auc**: 0.7550127730819666
- **depth_matched.per_peak_auc**: 0.6687866951547761
- **depth_matched.n_bins**: 5
- **depth_matched.cell_depth_median**: 14.0
