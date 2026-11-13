"""scatac-calib：scATAC 序列模型的校准评测命令行。

子命令:
  zero-model   零模型 AUROC（逐峰/逐细胞）
  oracle       群级 oracle（需要 cell -> group 标签）
  delta        ΔAUROC + bootstrap CI + 配对检验
  depth        深度分层后的逐峰 AUROC
  audit        一次跑完以上全部，并输出最低报告清单

零模型有两种口径（``--null``）：
  frequency    频率零模型。逐细胞轴用 test 峰检出率（transductive），
               逐峰轴用训练矩阵的细胞检出率；``--train-matrix`` 必填。
  inductive    5 折 cell cross-fitting。每折细胞的 null 分数只用其余折细胞的
               峰检出频率（被评测细胞的标签不参与），与工单 047 一致。
"""

import argparse
import json
import sys

import numpy as np

from . import baselines, bootstrap, io, metrics, report


def _fail(msg):
    """输入校验失败：写 stderr 并以退出码 2 结束（不是 warning）。"""
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(2)


def _common(p):
    p.add_argument("--matrix", required=True,
                   help="test 标签矩阵 (.npz/.h5ad，默认峰 x 细胞)")
    p.add_argument("--preds", default=None,
                   help="模型预测 (.npz 含 preds 键 / .npy)；zero-model/oracle 可省略")
    p.add_argument("--preds-key", default="preds")
    p.add_argument("--train-matrix", default=None,
                   help="训练矩阵 (峰 x 细胞)；frequency null 下必填")
    p.add_argument("--groups", default=None, help="test 细胞分组 (barcode<TAB>cluster)")
    p.add_argument("--train-groups", default=None,
                   help="训练细胞分组；缺省时复用 --groups")
    p.add_argument("--matrix-orientation", choices=io.ORIENTATIONS, default="auto",
                   help="矩阵方向：auto / peaks-x-cells / cells-x-peaks")
    p.add_argument("--null", choices=("frequency", "inductive"), default="frequency",
                   help="零模型口径；inductive 为 5-fold cell cross-fitting")
    p.add_argument("--n-folds", type=int, default=5,
                   help="inductive null 的折数，默认 5")
    p.add_argument("--fold-seed", type=int, default=20261112,
                   help="inductive null 的折划分种子，默认 20261112")
    p.add_argument("--out", required=True, help="输出目录")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261112)
    p.add_argument("--depth-bins", type=int, default=10)


def _load(args, need_preds=True):
    labels, minfo = io.load_matrix(
        args.matrix, orientation=args.matrix_orientation, return_info=True
    )
    labels = labels.toarray()
    n_peaks, n_cells = labels.shape

    preds = None
    if need_preds:
        if not args.preds:
            _fail("--preds is required for this command")
        preds = io.load_preds(args.preds, args.preds_key)
        if preds.shape != labels.shape:
            _fail(
                f"shape mismatch: preds {preds.shape} vs labels {labels.shape} "
                f"(matrix read as {minfo['orientation']}); if the matrix is "
                f"transposed, pass --matrix-orientation"
            )

    train = None
    if args.train_matrix:
        train, tinfo = io.load_matrix(
            args.train_matrix, orientation=args.matrix_orientation, return_info=True
        )
        if train.shape[1] != n_cells:
            _fail(
                f"--train-matrix has {train.shape[1]} cells but --matrix has "
                f"{n_cells}; the null needs matching cell columns "
                f"(train read as {tinfo['orientation']}; check --matrix-orientation)"
            )

    groups = io.load_groups(args.groups) if args.groups else None
    if groups is not None and len(groups) != n_cells:
        _fail(f"--groups has {len(groups)} entries but --matrix has {n_cells} cells")
    train_groups = io.load_groups(args.train_groups) if args.train_groups else groups
    if train is not None and train_groups is not None:
        if len(train_groups) != train.shape[1]:
            _fail(
                f"--train-groups has {len(train_groups)} entries but "
                f"--train-matrix has {train.shape[1]} cells"
            )
    return labels, preds, train, groups, train_groups, minfo


def _require_train(args, command):
    """frequency null 依赖训练矩阵；缺失是硬错误，不静默退化。"""
    if args.null == "frequency" and not args.train_matrix:
        _fail(
            f"{command}: --train-matrix is required when --null=frequency "
            f"(the frequency null needs training cell detection rates); "
            f"omit it only with --null=inductive"
        )


def _core_metrics(labels, preds):
    peak_auc = metrics.per_peak_auc(labels, preds)
    cell_auc = metrics.per_cell_auc(labels, preds)
    return {
        "metrics": metrics.mean_auc(peak_auc, cell_auc),
        "ap": metrics.average_precision(labels, preds),
        "_peak_auc": peak_auc,
        "_cell_auc": cell_auc,
    }


def _frequency_null(labels, train, require_train):
    """频率零模型；返回 (block, z_peak, z_cell)。

    ``require_train=True`` 时没有训练矩阵直接报错（CLI 的 frequency 口径）；
    inductive 模式下允许退化为 test 边际，但 ``cell_freq_source`` 会写明。
    """
    peak_freq = np.asarray(labels.mean(axis=1)).ravel().astype(np.float64)
    if train is not None:
        cell_freq = np.asarray(train.mean(axis=0)).ravel().astype(np.float64)
        source = "train"
    elif require_train:
        _fail("frequency null requires --train-matrix")
    else:
        cell_freq = np.asarray(labels.mean(axis=0)).ravel().astype(np.float64)
        source = "test(label-informed)"
    z_peak, z_cell = baselines.zero_auc_arrays(labels, peak_freq, cell_freq)
    block = {
        "zero_model": {
            "per_peak_auc": float(np.nanmean(z_peak)),
            "n_valid_peaks": int(np.sum(~np.isnan(z_peak))),
            "per_cell_auc": float(np.nanmean(z_cell)),
            "n_valid_cells": int(np.sum(~np.isnan(z_cell))),
            "cell_freq_source": source,
            "note": (
                "只用边际频率打分，不含序列信息；peak_freq 来自 test 峰检出率"
                "（逐细胞轴，transductive），cell_freq 来自 "
                + ("train 细胞检出率" if source == "train"
                   else "test 细胞检出率（label-informed，仅作参考）")
            ),
        },
        "_z_peak": z_peak,
        "_z_cell": z_cell,
    }
    return block, z_peak, z_cell


def _inductive_block(labels, core, groups, args, rng):
    """5-fold cell cross-fitting 的 inductive null 与 Δ。"""
    fold_ids, meta = baselines.make_cell_folds(
        labels.shape[1], groups=groups, n_folds=args.n_folds, seed=args.fold_seed
    )
    null_auc, per_fold = baselines.inductive_null_auc(labels, fold_ids)
    null_mean, n_null_valid = baselines.weighted_fold_mean(null_auc, fold_ids)
    block = {
        "inductive_null": {
            "kind": "5-fold cell cross-fitting",
            "per_cell_auc": float(null_mean),
            "n_valid_cells": int(n_null_valid),
            "n_folds": meta["n_folds"],
            "fold_seed": meta["fold_seed"],
            "fold_method": meta["method"],
            "stratified_by_groups": meta["stratified_by_groups"],
            "groups_missing_some_fold": meta["groups_missing_some_fold"],
            "fold_sizes": meta["fold_sizes"],
            "per_fold": per_fold,
            "null_score": (
                "peak detection frequency estimated from the other folds' cells; "
                "the evaluated cell's labels never enter its own null score"
            ),
            "note": (
                "逐细胞轴的 inductive 口径；逐峰轴仍用 train 的细胞检出率"
                "（对 test 峰已是 inductive）"
            ),
        }
    }
    if core is not None:
        model_auc = core["_cell_auc"]
        model_mean, n_model_valid = baselines.weighted_fold_mean(model_auc, fold_ids)
        paired = bootstrap.boot_paired_folds(
            model_auc, null_auc, fold_ids, rng, args.n_boot
        )
        block["delta_vs_inductive_null"] = {
            "cell": {
                "model": float(model_mean),
                "model_n_valid_cells": int(n_model_valid),
                "baseline": float(null_mean),
                "delta": float(model_mean - null_mean),
                "paired": paired,
            }
        }
    return block


def _oracle(labels, train, groups, train_groups):
    if groups is None:
        return {}
    # 峰集 split 时 train/test 的峰不是同一批，profile 必须来自同一批峰。
    # - 峰数一致：用 train（无泄漏，推荐）
    # - 峰数不一致：用 test 自身的群均值（label-informed 上界，需在文中注明）
    if train is not None and train.shape[0] == labels.shape[0] and train_groups is not None:
        src, src_groups, source = train, train_groups, "train"
    else:
        src, src_groups, source = labels, groups, "test(label-informed)"
    uniq, prof, index = baselines.group_profiles(src, src_groups, labels.shape[0])
    o_peak, o_cell, missing, _ = baselines.oracle_auc_arrays(
        labels, groups, uniq, prof, index
    )
    return {
        "oracle": {
            "kind": "group-level (cluster profile)",
            "source": source,
            "n_groups": len(uniq),
            "missing_cells": int(missing),
            "per_peak_auc": float(np.nanmean(o_peak)),
            "per_cell_auc": float(np.nanmean(o_cell)),
            "note": "已知细胞所属群时的上限；不代表模型可达到",
        },
        "_o_peak": o_peak,
        "_o_cell": o_cell,
    }


def _depth(labels, preds, args):
    depth = np.asarray(labels.sum(axis=0)).ravel()
    matched, per_bin = metrics.depth_matched_peak_auc(
        labels, preds, depth, args.depth_bins
    )
    return {
        "depth_matched": {
            "per_peak_auc": matched,
            "n_bins": args.depth_bins,
            "per_bin_auc": [float(x) for x in per_bin],
            "cell_depth_median": float(np.median(depth)),
        }
    }


def cmd_audit(args):
    _require_train(args, "audit")
    labels, preds, train, groups, train_groups, minfo = _load(args, need_preds=True)
    rng = np.random.default_rng(args.seed)
    out = io.ensure_dir(args.out)

    result = {
        "dataset": args.matrix,
        "preds": args.preds,
        "null": args.null,
        "matrix_format": minfo["format"],
        "matrix_orientation": minfo["orientation"],
        "n_peaks": int(labels.shape[0]),
        "n_cells": int(labels.shape[1]),
    }
    core = _core_metrics(labels, preds)
    result["metrics"] = core["metrics"]
    result["ap"] = core["ap"]

    freq, z_peak, z_cell = _frequency_null(labels, train, require_train=True)
    result["zero_model"] = freq["zero_model"]
    result["delta_vs_zero"] = {
        "peak": bootstrap.delta_summary(core["_peak_auc"], z_peak, rng, args.n_boot),
        "cell": bootstrap.delta_summary(core["_cell_auc"], z_cell, rng, args.n_boot),
    }
    result["ci95"] = {
        "per_peak_auc": bootstrap.boot_mean_ci(core["_peak_auc"], rng, args.n_boot),
        "per_cell_auc": bootstrap.boot_mean_ci(core["_cell_auc"], rng, args.n_boot),
    }

    if args.null == "inductive":
        # inductive 的 bootstrap 用独立的 RNG，避免改变 v1.0 frequency 口径的 CI。
        rng_ind = np.random.default_rng(args.seed)
        result.update(_inductive_block(labels, core, groups, args, rng_ind))

    o = _oracle(labels, train, groups, train_groups)
    if o:
        result["oracle"] = o["oracle"]

    result.update(_depth(labels, preds, args))

    report.write_json(out / "calib_report.json", _strip_private(result))
    report.write_checklist(out / "reporting_checklist.md", result, args.matrix)
    print(json.dumps(_strip_private(result), indent=2, ensure_ascii=False))
    return result


def cmd_zero(args):
    _require_train(args, "zero-model")
    labels, preds, train, groups, _, _ = _load(args, need_preds=False)
    out = io.ensure_dir(args.out)
    freq, _, _ = _frequency_null(
        labels, train, require_train=args.null == "frequency"
    )
    result = {"null": args.null, **freq}
    if args.null == "inductive":
        result.update(_inductive_block(labels, None, groups, args, None))
    report.write_json(out / "zero_model.json", _strip_private(result))
    if args.null == "inductive":
        report.write_json(out / "inductive_null.json",
                          _strip_private(result)["inductive_null"])
    print(json.dumps(_strip_private(result), indent=2, ensure_ascii=False))


def cmd_oracle(args):
    labels, preds, train, groups, train_groups, _ = _load(args, need_preds=False)
    out = io.ensure_dir(args.out)
    o = _oracle(labels, train, groups, train_groups)
    if not o:
        raise SystemExit("oracle 需要 --groups（以及可选的 --train-matrix）")
    report.write_json(out / "oracle.json", _strip_private(o))
    print(json.dumps(_strip_private(o), indent=2, ensure_ascii=False))


def cmd_delta(args):
    _require_train(args, "delta")
    labels, preds, train, groups, _, _ = _load(args, need_preds=True)
    out = io.ensure_dir(args.out)
    rng = np.random.default_rng(args.seed)
    core = _core_metrics(labels, preds)
    res = {}
    if train is not None:
        freq, z_peak, z_cell = _frequency_null(labels, train, require_train=False)
        res["zero_model"] = freq["zero_model"]
        res["delta_vs_zero"] = {
            "peak": bootstrap.delta_summary(core["_peak_auc"], z_peak, rng, args.n_boot),
            "cell": bootstrap.delta_summary(core["_cell_auc"], z_cell, rng, args.n_boot),
        }
    res["ci95"] = {
        "per_peak_auc": bootstrap.boot_mean_ci(core["_peak_auc"], rng, args.n_boot),
        "per_cell_auc": bootstrap.boot_mean_ci(core["_cell_auc"], rng, args.n_boot),
    }
    if args.null == "inductive":
        rng_ind = np.random.default_rng(args.seed)
        res.update(_inductive_block(labels, core, groups, args, rng_ind))
    report.write_json(out / "delta.json", _strip_private(res))
    print(json.dumps(_strip_private(res), indent=2, ensure_ascii=False))


def cmd_depth(args):
    labels, preds, _, _, _, _ = _load(args, need_preds=True)
    out = io.ensure_dir(args.out)
    res = _depth(labels, preds, args)
    report.write_json(out / "depth.json", res)
    print(json.dumps(res, indent=2, ensure_ascii=False))


def _strip_private(obj):
    return {k: v for k, v in obj.items() if not k.startswith("_")}


def build_parser():
    p = argparse.ArgumentParser(prog="scatac-calib", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name, fn, help_ in [
        ("audit", cmd_audit, "一次跑完零模型/oracle/delta/depth 并输出清单"),
        ("zero-model", cmd_zero, "零模型 AUROC（--null frequency/inductive）"),
        ("oracle", cmd_oracle, "群级 oracle（--preds 可省略）"),
        ("delta", cmd_delta, "ΔAUROC + bootstrap CI + 配对检验"),
        ("depth", cmd_depth, "深度分层逐峰 AUROC"),
    ]:
        sp = sub.add_parser(name, help=help_)
        _common(sp)
        sp.set_defaults(func=fn)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
