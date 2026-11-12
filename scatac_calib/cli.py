"""scatac-calib：scATAC 序列模型的校准评测命令行。

子命令:
  zero-model   零模型 AUROC（逐峰/逐细胞）
  oracle       群级 oracle（需要 cell -> group 标签）
  delta        ΔAUROC + bootstrap CI + 配对检验
  depth        深度分层后的逐峰 AUROC
  audit        一次跑完以上全部，并输出最低报告清单
"""

import argparse
import json
import sys

import numpy as np

from . import baselines, bootstrap, io, metrics, report


def _common(p):
    p.add_argument("--matrix", required=True, help="test 标签矩阵 (.npz/.h5ad，峰 x 细胞)")
    p.add_argument("--preds", required=True, help="模型预测 (.npz 含 preds 键 / .npy)")
    p.add_argument("--preds-key", default="preds")
    p.add_argument("--train-matrix", default=None, help="训练矩阵，用于零模型边际频率")
    p.add_argument("--groups", default=None, help="test 细胞分组 (barcode<TAB>cluster)")
    p.add_argument("--train-groups", default=None, help="训练细胞分组；缺省时复用 --groups")
    p.add_argument("--out", required=True, help="输出目录")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261112)
    p.add_argument("--depth-bins", type=int, default=10)


def _load(args):
    labels = io.load_matrix(args.matrix).toarray()
    preds = io.load_preds(args.preds, args.preds_key)
    if preds.shape != labels.shape:
        raise SystemExit(f"shape mismatch: preds {preds.shape} vs labels {labels.shape}")
    train = io.load_matrix(args.train_matrix) if args.train_matrix else None
    groups = io.load_groups(args.groups) if args.groups else None
    train_groups = io.load_groups(args.train_groups) if args.train_groups else groups
    return labels, preds, train, groups, train_groups


def _core_metrics(labels, preds):
    peak_auc = metrics.per_peak_auc(labels, preds)
    cell_auc = metrics.per_cell_auc(labels, preds)
    return {
        "metrics": metrics.mean_auc(peak_auc, cell_auc),
        "ap": metrics.average_precision(labels, preds),
        "_peak_auc": peak_auc,
        "_cell_auc": cell_auc,
    }


def _zero(labels, train):
    peak_freq, cell_freq, source = baselines.zero_scores(train, labels)
    z_peak, z_cell = baselines.zero_auc_arrays(labels, peak_freq, cell_freq)
    return {
        "zero_model": {
            "per_peak_auc": float(np.nanmean(z_peak)),
            "per_cell_auc": float(np.nanmean(z_cell)),
            "cell_freq_source": source,
            "note": "只用边际频率打分，不含序列信息；peak_freq 来自 test 峰检出率，cell_freq 来自 train 细胞检出率",
        },
        "_z_peak": z_peak,
        "_z_cell": z_cell,
    }


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
    labels, preds, train, groups, train_groups = _load(args)
    rng = np.random.default_rng(args.seed)
    out = io.ensure_dir(args.out)

    result = {"dataset": args.matrix, "preds": args.preds,
              "n_peaks": int(labels.shape[0]), "n_cells": int(labels.shape[1])}
    core = _core_metrics(labels, preds)
    result["metrics"] = core["metrics"]
    result["ap"] = core["ap"]

    z = _zero(labels, train)
    result["zero_model"] = z["zero_model"]
    result["delta_vs_zero"] = {
        "peak": bootstrap.delta_summary(core["_peak_auc"], z["_z_peak"], rng, args.n_boot),
        "cell": bootstrap.delta_summary(core["_cell_auc"], z["_z_cell"], rng, args.n_boot),
    }
    result["ci95"] = {
        "per_peak_auc": bootstrap.boot_mean_ci(core["_peak_auc"], rng, args.n_boot),
        "per_cell_auc": bootstrap.boot_mean_ci(core["_cell_auc"], rng, args.n_boot),
    }

    o = _oracle(labels, train, groups, train_groups)
    if o:
        result["oracle"] = o["oracle"]

    result.update(_depth(labels, preds, args))

    report.write_json(out / "calib_report.json", _strip_private(result))
    report.write_checklist(out / "reporting_checklist.md", result, args.matrix)
    print(json.dumps(_strip_private(result), indent=2, ensure_ascii=False))
    return result


def cmd_zero(args):
    labels, preds, train, _, _ = _load(args)
    out = io.ensure_dir(args.out)
    z = _zero(labels, train)
    report.write_json(out / "zero_model.json", _strip_private(z))
    print(json.dumps(_strip_private(z), indent=2, ensure_ascii=False))


def cmd_oracle(args):
    labels, preds, train, groups, train_groups = _load(args)
    out = io.ensure_dir(args.out)
    o = _oracle(labels, train, groups, train_groups)
    if not o:
        raise SystemExit("oracle 需要 --train-matrix、--groups（以及匹配的 train 细胞）")
    report.write_json(out / "oracle.json", _strip_private(o))
    print(json.dumps(_strip_private(o), indent=2, ensure_ascii=False))


def cmd_delta(args):
    labels, preds, train, _, _ = _load(args)
    out = io.ensure_dir(args.out)
    rng = np.random.default_rng(args.seed)
    core = _core_metrics(labels, preds)
    z = _zero(labels, train)
    res = {
        "delta_vs_zero": {
            "peak": bootstrap.delta_summary(core["_peak_auc"], z["_z_peak"], rng, args.n_boot),
            "cell": bootstrap.delta_summary(core["_cell_auc"], z["_z_cell"], rng, args.n_boot),
        },
        "ci95": {
            "per_peak_auc": bootstrap.boot_mean_ci(core["_peak_auc"], rng, args.n_boot),
            "per_cell_auc": bootstrap.boot_mean_ci(core["_cell_auc"], rng, args.n_boot),
        },
    }
    report.write_json(out / "delta.json", res)
    print(json.dumps(res, indent=2, ensure_ascii=False))


def cmd_depth(args):
    labels, preds, _, _, _ = _load(args)
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
        ("zero-model", cmd_zero, "零模型 AUROC"),
        ("oracle", cmd_oracle, "群级 oracle"),
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
