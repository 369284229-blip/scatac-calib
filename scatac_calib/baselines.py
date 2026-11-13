"""零模型与群级 oracle。

零模型（frequency null）：只用训练集的边际频率给 test 打分，不含任何序列信息。
  - cell_freq: 每个峰在训练细胞中的可及比例 -> 用来给“峰 i 的细胞轴”打分
  - peak_freq: 每个细胞在训练峰中的可及比例 -> 用来给“细胞 j 的峰轴”打分

群级 oracle（group-level oracle）：假设已知每个细胞所属的群（如 Leiden cluster），
用该群的峰可及性均值作打分。它代表“只知道细胞类型、不知道单细胞”时的上限。
"""

import numpy as np
from sklearn.model_selection import KFold, StratifiedKFold

from .metrics import safe_auc


def marginal_frequencies(train_matrix):
    """train_matrix: CSR (n_peaks, n_cells) -> peak_freq(n_cells), cell_freq(n_peaks)"""
    peak_freq = np.asarray(train_matrix.mean(axis=0)).ravel().astype(np.float64)
    cell_freq = np.asarray(train_matrix.mean(axis=1)).ravel().astype(np.float64)
    return peak_freq, cell_freq


def zero_scores(train_matrix, labels):
    """返回对 test 打分的两个边际频率（口径与 minloop_eval_wheat.py 一致）。

    - peak_freq: 每个 test 峰在 test 细胞中的检出比例（峰轴打分，逐细胞 AUC）
    - cell_freq: 每个细胞在 train 峰中的检出比例（细胞轴打分，逐峰 AUC）
      没有 train 矩阵或细胞数不匹配时退化为 test 边际。

    注意：峰集 split 时 train/test 峰不同，cell_freq 只能来自 train 的列均值，
    这正是项目里 S-013/S-014 的口径。
    """
    n_peaks, n_cells = labels.shape
    peak_freq = np.asarray(labels.mean(axis=1)).ravel().astype(np.float64)
    if train_matrix is not None and train_matrix.shape[1] == n_cells:
        cell_freq = np.asarray(train_matrix.mean(axis=0)).ravel().astype(np.float64)
        source = "train"
    else:
        cell_freq = np.asarray(labels.mean(axis=0)).ravel().astype(np.float64)
        source = "test"
    return peak_freq, cell_freq, source


def zero_auc_arrays(labels, peak_freq, cell_freq):
    """逐峰轴用 cell_freq（每个细胞一个分数）打分；逐细胞轴用 peak_freq（每个峰一个分数）打分。"""
    n_peaks = labels.shape[0]
    z_peak = np.full(n_peaks, np.nan)
    for i in range(n_peaks):
        z_peak[i] = safe_auc(labels[i], cell_freq)
    z_cell = np.full(labels.shape[1], np.nan)
    for j in range(labels.shape[1]):
        z_cell[j] = safe_auc(labels[:, j], peak_freq)
    return z_peak, z_cell


def make_cell_folds(n_cells, groups=None, n_folds=5, seed=20261112):
    """构造 cell cross-fitting 的折。

    ``groups`` 给定时按群分层（与工单 047 的
    ``StratifiedKFold(5, shuffle=True, random_state=20261112)`` 一致）；
    否则退化为不分层 KFold，并在 ``method`` 里注明。

    返回 ``(fold_ids, meta)``：``fold_ids[k]`` 是第 k 折的细胞列索引。
    """
    if n_folds < 2:
        raise ValueError("n_folds must be >= 2")
    folds = min(int(n_folds), int(n_cells))
    if folds < 2:
        raise ValueError(
            f"cell cross-fitting needs at least 2 cells, got {n_cells}"
        )
    if groups is not None:
        groups = np.asarray(groups)
        if len(groups) != n_cells:
            raise ValueError(
                f"groups has {len(groups)} entries but matrix has {n_cells} cells"
            )
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
        fold_ids = [te for _, te in splitter.split(np.zeros(n_cells), groups)]
        method = "stratified-by-groups"
        unique = np.unique(groups)
        absent = {
            str(g): [k for k, idx in enumerate(fold_ids)
                     if g not in set(groups[idx].tolist())]
            for g in unique
        }
        groups_missing = {g: ks for g, ks in absent.items() if ks}
    else:
        splitter = KFold(n_splits=folds, shuffle=True, random_state=seed)
        fold_ids = [te for _, te in splitter.split(np.zeros(n_cells))]
        method = "kfold-unstratified"
        groups_missing = {}
    meta = {
        "n_folds": folds,
        "fold_seed": int(seed),
        "method": method,
        "stratified_by_groups": groups is not None,
        "groups_missing_some_fold": groups_missing,
        "fold_sizes": [int(len(idx)) for idx in fold_ids],
    }
    return fold_ids, meta


def inductive_null_auc(labels, fold_ids):
    """逐细胞 inductive null：每折用其余折的峰检出频率给本折细胞打分。

    被评测细胞的标签不进入自己的 null 分数。返回 ``(null_auc, per_fold)``；
    ``null_auc[j]`` 是细胞 j 的 null AUC（单类细胞为 NaN）。
    """
    n_cells = labels.shape[1]
    null_auc = np.full(n_cells, np.nan, dtype=np.float64)
    per_fold = []
    for k, test_idx in enumerate(fold_ids):
        train_mask = np.ones(n_cells, dtype=bool)
        train_mask[test_idx] = False
        peak_freq = np.asarray(labels[:, train_mask].mean(axis=1)).ravel()
        peak_freq = peak_freq.astype(np.float64)
        for j in test_idx:
            null_auc[j] = safe_auc(labels[:, j], peak_freq)
        block = null_auc[test_idx]
        per_fold.append({
            "fold": int(k),
            "n_cells": int(len(test_idx)),
            "n_valid_cells": int(np.sum(~np.isnan(block))),
            "null_cell_auc": float(np.nanmean(block)),
        })
    return null_auc, per_fold


def weighted_fold_mean(values, fold_ids):
    """按折内有效细胞数加权汇总（与工单 047 的 ``weighted_mean`` 一致）。"""
    values = np.asarray(values, dtype=np.float64)
    total = 0.0
    n = 0
    for idx in fold_ids:
        v = values[idx]
        v = v[~np.isnan(v)]
        total += float(v.sum())
        n += len(v)
    return (total / n if n else float("nan")), n


def group_profiles(train_matrix, train_groups, n_peaks):
    """群 x 峰 的可及性均值。train_matrix: (n_peaks_train, n_cells_train)"""
    uniq = sorted(set(train_groups.tolist()))
    prof = np.zeros((len(uniq), n_peaks), dtype=np.float32)
    index = {g: k for k, g in enumerate(uniq)}
    cols_by_group = {g: [] for g in uniq}
    for j, g in enumerate(train_groups):
        cols_by_group[g].append(j)
    for g, cols in cols_by_group.items():
        prof[index[g]] = np.asarray(train_matrix[:, cols].mean(axis=1)).ravel()
    return uniq, prof, index


def oracle_auc_arrays(test_labels, test_groups, uniq, prof, index):
    """群级 oracle：test 细胞按自己群的 profile 打分。"""
    n_peaks, n_cells = test_labels.shape
    score = np.zeros((n_peaks, n_cells), dtype=np.float32)
    missing = 0
    for j, g in enumerate(test_groups):
        if g not in index:
            missing += 1
            score[:, j] = nan_profile(prof)
        else:
            score[:, j] = prof[index[g]]
    o_peak = np.full(n_peaks, np.nan)
    for i in range(n_peaks):
        o_peak[i] = safe_auc(test_labels[i], score[i])
    o_cell = np.full(n_cells, np.nan)
    for j in range(n_cells):
        o_cell[j] = safe_auc(test_labels[:, j], score[:, j])
    return o_peak, o_cell, missing, score


def nan_profile(prof):
    # 未知群退化为全体均值
    return prof.mean(axis=0)
