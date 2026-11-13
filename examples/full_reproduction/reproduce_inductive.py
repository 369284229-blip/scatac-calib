#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""在全量小麦 test split 上复现工单 047 的 inductive null 与 Δcell。

这不是打包进来的数据 fixture：完整 test 标签矩阵 4,665 × 8,661 与三个 seed
的预测共约 430 MB，不适合放进 wheel/sdist。本脚本只读项目数据，按 CLI 的同一
实现（``scatac_calib.baselines``）计算并核对：

    null（5-fold cell cross-fitting）  = 0.6852978   （容差 1e-3）
    Δcell（peak-selected 3 seed 均值） = +0.0193566 （容差 1e-3）

用法（数据在服务器 L-wdl 上时）：

    python examples/full_reproduction/reproduce_inductive.py \
      --data-dir /root/wheat_project/02_results/wheat/dataset_v3 \
      --preds \
        /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed0/preds.npz \
        /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed1/preds.npz \
        /root/wheat_project/02_results/wheat/p1_benchmark/bootstrap/seed2/preds.npz

退出码：0 = 两个目标都通过；1 = 有任一项超出容差。
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# 允许从仓库根目录直接运行
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scatac_calib import baselines, io, metrics  # noqa: E402

EXPECTED_NULL = 0.6852978371234296
EXPECTED_DELTA_MEAN = 0.0193566


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", required=True,
                   help="含 m_test.npz / cell_groups.txt 的 dataset_v3 目录")
    p.add_argument("--matrix", default=None,
                   help="默认 <data-dir>/m_test.npz")
    p.add_argument("--groups", default=None,
                   help="默认 <data-dir>/cell_groups.txt")
    p.add_argument("--preds", nargs="+", required=True,
                   help="一个或多个 seed 的 preds.npz（3 seed 均值=工单 047 口径）")
    p.add_argument("--preds-key", default="preds")
    p.add_argument("--n-folds", type=int, default=5)
    p.add_argument("--fold-seed", type=int, default=20261112)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--out", default=None, help="可选：写出 JSON 结果")
    args = p.parse_args(argv)

    data = Path(args.data_dir)
    matrix_path = Path(args.matrix) if args.matrix else data / "m_test.npz"
    groups_path = Path(args.groups) if args.groups else data / "cell_groups.txt"
    labels = io.load_matrix(matrix_path).toarray()
    groups = io.load_groups(groups_path)
    if len(groups) != labels.shape[1]:
        raise SystemExit(
            f"groups {len(groups)} != cells {labels.shape[1]}"
        )

    folds, meta = baselines.make_cell_folds(
        labels.shape[1], groups=groups, n_folds=args.n_folds, seed=args.fold_seed
    )
    null_auc, per_fold = baselines.inductive_null_auc(labels, folds)
    null_mean, n_valid = baselines.weighted_fold_mean(null_auc, folds)

    per_seed = []
    for path in args.preds:
        preds = io.load_preds(path, args.preds_key)
        if preds.shape != labels.shape:
            raise SystemExit(
                f"{path}: preds shape {preds.shape} != labels {labels.shape}"
            )
        model_auc = metrics.per_cell_auc(labels, preds)
        model_mean, n_model = baselines.weighted_fold_mean(model_auc, folds)
        per_seed.append({
            "preds": str(path),
            "model_cell_auc": float(model_mean),
            "n_valid_cells": int(n_model),
            "delta": float(model_mean - null_mean),
        })
    delta_mean = float(np.mean([x["delta"] for x in per_seed]))

    result = {
        "n_peaks": int(labels.shape[0]),
        "n_cells": int(labels.shape[1]),
        "folds": meta,
        "null_cell_auc_inductive": float(null_mean),
        "n_valid_cells": int(n_valid),
        "per_fold": per_fold,
        "per_seed": per_seed,
        "delta_mean": delta_mean,
        "expected": {
            "null_cell_auc_inductive": EXPECTED_NULL,
            "delta_mean": EXPECTED_DELTA_MEAN,
            "tolerance": args.tol,
        },
        "pass": {
            "null": bool(abs(null_mean - EXPECTED_NULL) <= args.tol),
            "delta_mean": bool(abs(delta_mean - EXPECTED_DELTA_MEAN) <= args.tol),
        },
    }

    print(json.dumps(result, indent=2, ensure_ascii=False))
    if args.out:
        io.ensure_dir(Path(args.out).parent)
        Path(args.out).write_text(json.dumps(result, indent=2, ensure_ascii=False))
    if not all(result["pass"].values()):
        print("FAIL: value(s) outside tolerance", file=sys.stderr)
        return 1
    print("PASS: inductive null and delta reproduce work order 047", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
