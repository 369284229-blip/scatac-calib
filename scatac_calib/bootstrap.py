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
