"""逐峰 / 逐细胞 AUROC、AP，以及深度分层控制。"""

import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

EPS = 1e-12


def safe_auc(y, s):
    """单一类标签返回 NaN（不参与平均）。"""
    if y.min() == y.max():
        return np.nan
    return float(roc_auc_score(y, s))


def per_peak_auc(labels, preds, idx=None):
    """每一行（峰）跨细胞的 AUROC。"""
    n = labels.shape[0] if idx is None else len(idx)
    out = np.full(n, np.nan, dtype=np.float64)
    rows = range(labels.shape[0]) if idx is None else idx
    for k, i in enumerate(rows):
        out[k] = safe_auc(labels[i], preds[i])
    return out


def per_cell_auc(labels, preds, idx=None):
    """每一列（细胞）跨峰的 AUROC。"""
    out = np.full(labels.shape[1], np.nan, dtype=np.float64)
    cols = range(labels.shape[1]) if idx is None else idx
    for j in cols:
        out[j] = safe_auc(labels[:, j], preds[:, j])
    return out


def mean_auc(peak_auc=None, cell_auc=None):
    res = {}
    if peak_auc is not None:
        res["per_peak_auc"] = float(np.nanmean(peak_auc))
        res["n_valid_peaks"] = int(np.sum(~np.isnan(peak_auc)))
    if cell_auc is not None:
        res["per_cell_auc"] = float(np.nanmean(cell_auc))
        res["n_valid_cells"] = int(np.sum(~np.isnan(cell_auc)))
    return res


def average_precision(labels, preds):
    return float(average_precision_score(labels.flatten(), preds.flatten()))


def depth_bins(depth, n_bins=10):
    """按细胞深度分位数返回每个细胞所属的 bin。"""
    depth = np.asarray(depth, dtype=np.float64)
    # rank 后等分，避免并列值把某些 bin 抽空
    order = np.argsort(depth, kind="mergesort")
    ranks = np.empty_like(order)
    ranks[order] = np.arange(len(depth))
    return np.minimum((ranks * n_bins) // max(len(depth), 1), n_bins - 1)


def depth_matched_peak_auc(labels, preds, depth, n_bins=10):
    """在每个深度 bin 内算逐峰 AUC，再对峰取平均。

    回答：模型的逐峰表现有多少来自“深细胞更容易被判对”。
    """
    bins = depth_bins(depth, n_bins)
    n_peaks = labels.shape[0]
    per_bin = []
    for b in np.unique(bins):
        cells = np.where(bins == b)[0]
        if len(cells) < 2:
            continue
        aucs = np.full(n_peaks, np.nan)
        for i in range(n_peaks):
            y = labels[i, cells]
            if y.min() == y.max():
                continue
            aucs[i] = roc_auc_score(y, preds[i, cells])
        if not np.all(np.isnan(aucs)):
            per_bin.append(np.nanmean(aucs))
    return float(np.mean(per_bin)) if per_bin else float("nan"), per_bin
