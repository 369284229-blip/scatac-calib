"""Inductive（5-fold cell cross-fitting）零模型的单元与端到端测试。"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scatac_calib import baselines  # noqa: E402

FIXTURE = ROOT / "examples" / "real_fixture"


def _run_cli(args, cwd):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "scatac_calib", *args],
        cwd=cwd, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )


class TestFoldSplit(unittest.TestCase):
    def test_folds_partition_all_cells_and_stratify(self):
        rng = np.random.default_rng(0)
        groups = np.repeat(np.arange(4), 25)
        fold_ids, meta = baselines.make_cell_folds(
            len(groups), groups=groups, n_folds=5, seed=20261112
        )
        all_idx = np.sort(np.concatenate(fold_ids))
        np.testing.assert_array_equal(all_idx, np.arange(len(groups)))
        self.assertEqual(len(fold_ids), 5)
        self.assertEqual(meta["method"], "stratified-by-groups")
        self.assertTrue(meta["stratified_by_groups"])
        for idx in fold_ids:
            self.assertEqual(set(groups[idx].tolist()), {0, 1, 2, 3})

    def test_unstratified_fallback(self):
        fold_ids, meta = baselines.make_cell_folds(
            30, groups=None, n_folds=5, seed=1
        )
        self.assertEqual(meta["method"], "kfold-unstratified")
        self.assertEqual(sorted(np.concatenate(fold_ids).tolist()), list(range(30)))


class TestCrossFitting(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(7)
        self.labels = (rng.random((60, 40)) < 0.15).astype(np.uint8)
        self.groups = np.repeat(np.arange(4), 10)
        self.fold_ids, _ = baselines.make_cell_folds(
            40, groups=self.groups, n_folds=5, seed=20261112
        )

    def test_evaluated_fold_labels_do_not_enter_null_scores(self):
        """折内细胞的标签不应进入本折的 null 分数（关键科学性质）。

        直接截获传给 ``safe_auc`` 的 score 向量：改变折 0 的标签后，
        折 0 细胞拿到的 score（峰检出频率）必须逐位不变。
        """
        fold0 = self.fold_ids[0]
        calls = []
        original = baselines.safe_auc

        def spy(y, s):
            calls.append((np.array(y), np.array(s)))
            return original(y, s)

        baselines.safe_auc = spy
        try:
            baselines.inductive_null_auc(self.labels, self.fold_ids)
            base_calls = len(calls)
            self.assertEqual(base_calls, self.labels.shape[1])
            base_scores = [calls[j][1] for j in range(len(fold0))]
            base_y = [calls[j][0] for j in range(len(fold0))]
            base_score_fold1 = calls[len(fold0)][1].copy()

            flipped = self.labels.copy()
            flipped[:, fold0] = 1 - flipped[:, fold0]
            calls.clear()
            baselines.inductive_null_auc(flipped, self.fold_ids)
            new_scores = [calls[j][1] for j in range(len(fold0))]
            new_y = [calls[j][0] for j in range(len(fold0))]
        finally:
            baselines.safe_auc = original

        for a, b in zip(base_scores, new_scores):
            np.testing.assert_array_equal(a, b)
        # y 确实变了，说明上面的“不变”不是因为输入没变
        self.assertTrue(any(not np.array_equal(a, b) for a, b in zip(base_y, new_y)))
        # 其他折把折 0 当训练细胞，score 应该改变
        self.assertFalse(
            np.array_equal(calls[len(fold0)][1], base_score_fold1)
        )

    def test_null_scores_match_direct_computation(self):
        null_auc, per_fold = baselines.inductive_null_auc(
            self.labels, self.fold_ids
        )
        fold0 = self.fold_ids[0]
        mask = np.ones(self.labels.shape[1], dtype=bool)
        mask[fold0] = False
        freq = self.labels[:, mask].mean(axis=1)
        for j in fold0[:3]:
            y = self.labels[:, j]
            if y.min() == y.max():
                self.assertTrue(np.isnan(null_auc[j]))
                continue
            from sklearn.metrics import roc_auc_score
            self.assertAlmostEqual(
                null_auc[j], roc_auc_score(y, freq), delta=1e-12
            )
        self.assertEqual(len(per_fold), len(self.fold_ids))
        self.assertEqual(sum(x["n_cells"] for x in per_fold), 40)

    def test_weighted_fold_mean_matches_pooled_mean(self):
        values = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        folds = [np.array([0, 1]), np.array([2, 3, 4]), np.array([5, 6, 7, 8, 9])]
        mean, n = baselines.weighted_fold_mean(values, folds)
        self.assertAlmostEqual(mean, np.nanmean(values), delta=1e-12)
        self.assertEqual(n, 9)


class TestInductiveAuditFixture(unittest.TestCase):
    """在真实数据子集上端到端回归（方法口径，不是论文全量数字）。"""

    def test_audit_inductive_point_estimates(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            proc = _run_cli([
                "audit",
                "--matrix", str(FIXTURE / "m_test_fixture.npz"),
                "--preds", str(FIXTURE / "preds_fixture.npz"),
                "--train-matrix", str(FIXTURE / "m_train_fixture.npz"),
                "--groups", str(FIXTURE / "groups_fixture.txt"),
                "--null", "inductive",
                "--out", str(out),
                "--n-boot", "50",
                "--depth-bins", "5",
                "--seed", "20261112",
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            r = json.loads((out / "calib_report.json").read_text())
            self.assertEqual(r["null"], "inductive")
            self.assertEqual(r["inductive_null"]["n_folds"], 5)
            self.assertEqual(r["inductive_null"]["fold_seed"], 20261112)
            self.assertEqual(r["inductive_null"]["fold_sizes"], [16, 16, 16, 16, 16])
            self.assertEqual(r["inductive_null"]["n_valid_cells"], 64)
            self.assertAlmostEqual(
                r["inductive_null"]["per_cell_auc"], 0.6676299149688951, delta=1e-6
            )
            delta = r["delta_vs_inductive_null"]["cell"]
            self.assertAlmostEqual(delta["model"], 0.6940874244675193, delta=1e-6)
            self.assertAlmostEqual(delta["delta"], 0.026457509498624132, delta=1e-6)
            # 逐峰轴不受 null 口径影响
            self.assertAlmostEqual(
                r["zero_model"]["per_peak_auc"], 0.6834153421349665, delta=1e-6
            )
            self.assertEqual(r["zero_model"]["cell_freq_source"], "train")


if __name__ == "__main__":
    unittest.main()
