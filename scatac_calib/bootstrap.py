"""Bootstrap CI 与配对检验（口径与 项目1/代码/bootstrap_ci.py 一致）。"""

import numpy as np


def boot_mean_ci(x, rng, n_boot=1000):
    x = np.asarray(x, dtype=np.float64)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return [float("nan"), float("nan")]
    vals = np.empty(n_boot)
    for b in range(n_boot):
        vals[b] = x[rng.integers(0, len(x), len(x))].mean()
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def boot_paired(a, b, rng, n_boot=1000):
    """配对 bootstrap：mean(a-b) 的 95% CI 与双侧 p。"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ok = ~np.isnan(a) & ~np.isnan(b)
    a, b = a[ok], b[ok]
    if len(a) == 0:
        return {"ci95": [float("nan"), float("nan")], "p": float("nan"),
                "mean": float("nan"), "n": 0}
    deltas = np.empty(n_boot)
    for k in range(n_boot):
        idx = rng.integers(0, len(a), len(a))
        deltas[k] = a[idx].mean() - b[idx].mean()
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    p = 2.0 * min((deltas <= 0).mean(), (deltas >= 0).mean())
    return {"ci95": [float(lo), float(hi)], "p": float(min(1.0, p)),
            "mean": float(deltas.mean()), "n": int(len(a)), "n_boot": n_boot}


def boot_paired_folds(a, b, fold_ids, rng, n_boot=1000):
    """折内配对 bootstrap（口径与工单 047 一致）。

    每折有放回重采样同样数量的细胞，保留折结构；按折细胞数加权汇总
    ``mean(a-b)``。返回 95% CI、双侧 p 与均值。
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    diffs = a - b
    per_fold = []
    n_total = 0
    for idx in fold_ids:
        d = diffs[idx]
        d = d[~np.isnan(d)]
        per_fold.append((idx, d))
        n_total += len(d)
    if n_total == 0:
        return {"ci95": [float("nan"), float("nan")], "p": float("nan"),
                "mean": float("nan"), "n": 0, "n_boot": n_boot,
                "scheme": "per-fold cell bootstrap"}
    vals = np.empty(n_boot, dtype=np.float64)
    for k in range(n_boot):
        s = 0.0
        for idx, d in per_fold:
            if len(d) == 0:
                continue
            # 在有效细胞对应的折内位置重采样
            samp = rng.integers(0, len(d), len(d))
            s += len(d) * float(d[samp].mean())
        vals[k] = s / n_total
    lo, hi = np.percentile(vals, [2.5, 97.5])
    p = 2.0 * min((vals <= 0).mean(), (vals >= 0).mean())
    return {"ci95": [float(lo), float(hi)], "p": float(min(1.0, p)),
            "mean": float(vals.mean()), "n": int(n_total), "n_boot": n_boot,
            "scheme": "per-fold cell bootstrap"}


def delta_summary(a, b, rng, n_boot=1000):
    """a=模型, b=零模型/oracle；返回均值、CI、配对检验。"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ok = ~np.isnan(a) & ~np.isnan(b)
    return {
        "model": float(np.nanmean(a[ok])),
        "baseline": float(np.nanmean(b[ok])),
        "delta": float(np.nanmean(a[ok]) - np.nanmean(b[ok])),
        "paired": boot_paired(a, b, rng, n_boot),
    }
