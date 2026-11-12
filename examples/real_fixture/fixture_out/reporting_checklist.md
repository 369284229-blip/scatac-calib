# 最低报告清单：m_test_fixture.npz

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
- **metrics.per_peak_auc**: 0.7478242362161837
- **metrics.n_valid_peaks**: 73
- **metrics.per_cell_auc**: 0.6940874244675193
- **metrics.n_valid_cells**: 64
- **ap**: 0.07747587158056138
- **zero_model.per_peak_auc**: 0.6834153421349665
- **zero_model.per_cell_auc**: 0.8089406162485495
- **delta_vs_zero.peak.model**: 0.7478242362161837
- **delta_vs_zero.peak.baseline**: 0.6834153421349665
- **delta_vs_zero.peak.delta**: 0.0644088940812172
- **delta_vs_zero.peak.paired.p**: 0.0
- **delta_vs_zero.peak.paired.mean**: 0.06396887576445234
- **delta_vs_zero.peak.paired.n**: 73
- **delta_vs_zero.peak.paired.n_boot**: 1000
- **delta_vs_zero.cell.model**: 0.6940874244675193
- **delta_vs_zero.cell.baseline**: 0.8089406162485495
- **delta_vs_zero.cell.delta**: -0.11485319178103026
- **delta_vs_zero.cell.paired.p**: 0.0
- **delta_vs_zero.cell.paired.mean**: -0.11470921933719674
- **delta_vs_zero.cell.paired.n**: 64
- **delta_vs_zero.cell.paired.n_boot**: 1000
- **oracle.n_groups**: 14
- **oracle.missing_cells**: 0
- **oracle.per_peak_auc**: 0.9348451780725647
- **oracle.per_cell_auc**: 0.9550075335537698
- **depth_matched.per_peak_auc**: 0.6235446510334903
- **depth_matched.n_bins**: 5
- **depth_matched.cell_depth_median**: 2.0
