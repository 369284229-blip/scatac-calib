#!/usr/bin/env python3
"""Generate the synthetic demo inputs for scatac-calib.

The data are deliberately small and fully synthetic: 100 peaks x 80 cells,
4 cell groups with 25 group-marker peaks each.  They contain no real
sequence or expression data, so the demo can be shipped with the package.

Usage:
    python make_demo_data.py [--outdir DIR] [--seed INT]

Outputs (all written to ``--outdir``, default: this directory):
    m_demo.npz         test matrix, CSR, 100 peaks x 80 cells
    m_train_demo.npz   independent train matrix, same shape
    preds_demo.npz     model-like prediction matrix, key ``preds``
    groups_demo.txt    barcode<TAB>cluster, 80 rows
"""

import argparse
from pathlib import Path

import numpy as np
from scipy import sparse


def _rate_matrix(rng, n_peaks, n_cells, n_groups):
    cells_per_group = n_cells // n_groups
    base = rng.uniform(0.04, 0.12, size=(n_peaks, 1))
    group = np.repeat(np.arange(n_groups), cells_per_group)
    if len(group) < n_cells:  # tolerate cells not divisible by groups
        group = np.concatenate([group, np.zeros(n_cells - len(group), dtype=int)])

    rate = np.repeat(base, n_cells, axis=1).astype(np.float64)
    marker_idx = np.array_split(rng.permutation(n_peaks)[: n_peaks // 2], n_groups)
    for g, peaks in enumerate(marker_idx):
        rate[np.ix_(peaks, np.where(group == g)[0])] += 0.55
    return np.clip(rate, 0.0, 0.97), group


def _csr_npz(path, matrix):
    m = sparse.csr_matrix(matrix)
    np.savez_compressed(
        path,
        data=m.data,
        indices=m.indices,
        indptr=m.indptr,
        shape=np.asarray(m.shape, dtype=np.int64),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--seed", type=int, default=20261112)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    n_peaks, n_cells, n_groups = 100, 80, 4
    rate, group = _rate_matrix(rng, n_peaks, n_cells, n_groups)

    labels = (rng.random((n_peaks, n_cells)) < rate).astype(np.float32)
    train = (rng.random((n_peaks, n_cells)) < rate).astype(np.float32)

    # Model-like predictions: true rate plus independent noise.
    preds = np.clip(rate + rng.normal(0.0, 0.10, size=rate.shape), 1e-4, 1 - 1e-4)
    preds = preds.astype(np.float32)

    _csr_npz(out / "m_demo.npz", labels)
    _csr_npz(out / "m_train_demo.npz", train)
    np.savez_compressed(out / "preds_demo.npz", preds=preds)

    names = [f"G{g + 1}" for g in range(n_groups)]
    with open(out / "groups_demo.txt", "w") as f:
        for j in range(n_cells):
            f.write(f"cell_{j:03d}\t{names[group[j]]}\n")

    print(f"wrote {n_peaks} peaks x {n_cells} cells to {out}")
    print(f"label rate: {labels.mean():.4f}")
    print(f"groups: {dict(zip(*np.unique([names[g] for g in group], return_counts=True)))}")


if __name__ == "__main__":
    main()
