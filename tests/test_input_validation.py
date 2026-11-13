"""输入硬校验：--train-matrix、.h5ad 方向、--preds optional。"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scatac_calib import io  # noqa: E402

FIXTURE = ROOT / "examples" / "real_fixture"
DEMO = ROOT / "examples" / "demo"


def _run_cli(args, cwd):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "scatac_calib", *args],
        cwd=cwd, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )


def _save_csr(path, matrix):
    m = sparse.csr_matrix(matrix)
    np.savez_compressed(
        path, data=m.data, indices=m.indices, indptr=m.indptr, shape=m.shape
    )


class TestTrainMatrixRequired(unittest.TestCase):
    def test_frequency_null_without_train_matrix_hard_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = _run_cli([
                "audit",
                "--matrix", str(FIXTURE / "m_test_fixture.npz"),
                "--preds", str(FIXTURE / "preds_fixture.npz"),
                "--groups", str(FIXTURE / "groups_fixture.txt"),
                "--out", str(Path(tmp) / "out"),
                "--n-boot", "10",
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 2, msg=proc.stderr)
            self.assertIn("--train-matrix is required", proc.stderr)
            self.assertFalse((Path(tmp) / "out").exists())

    def test_frequency_zero_model_without_train_matrix_hard_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = _run_cli([
                "zero-model",
                "--matrix", str(FIXTURE / "m_test_fixture.npz"),
                "--out", str(Path(tmp) / "out"),
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--train-matrix is required", proc.stderr)

    def test_train_matrix_cell_count_mismatch_hard_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bad = tmp / "train_bad.npz"
            _save_csr(bad, sparse.random(
                20, 7, density=0.2, format="csr", random_state=0
            ))
            proc = _run_cli([
                "zero-model",
                "--matrix", str(FIXTURE / "m_test_fixture.npz"),
                "--train-matrix", str(bad),
                "--out", str(tmp / "out"),
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--train-matrix has 7 cells", proc.stderr)

    def test_inductive_allows_missing_train_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = _run_cli([
                "zero-model",
                "--matrix", str(FIXTURE / "m_test_fixture.npz"),
                "--groups", str(FIXTURE / "groups_fixture.txt"),
                "--null", "inductive",
                "--out", str(Path(tmp) / "out"),
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            r = json.loads((Path(tmp) / "out" / "zero_model.json").read_text())
            self.assertEqual(r["inductive_null"]["n_folds"], 5)
            self.assertEqual(r["zero_model"]["cell_freq_source"],
                             "test(label-informed)")


class TestPredsOptional(unittest.TestCase):
    def test_zero_model_runs_without_preds(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = _run_cli([
                "zero-model",
                "--matrix", str(DEMO / "m_demo.npz"),
                "--train-matrix", str(DEMO / "m_train_demo.npz"),
                "--out", str(Path(tmp) / "out"),
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            r = json.loads((Path(tmp) / "out" / "zero_model.json").read_text())
            self.assertIn("zero_model", r)
            self.assertNotIn("metrics", r)

    def test_oracle_runs_without_preds(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = _run_cli([
                "oracle",
                "--matrix", str(DEMO / "m_demo.npz"),
                "--train-matrix", str(DEMO / "m_train_demo.npz"),
                "--groups", str(DEMO / "groups_demo.txt"),
                "--out", str(Path(tmp) / "out"),
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            r = json.loads((Path(tmp) / "out" / "oracle.json").read_text())
            self.assertEqual(r["oracle"]["n_groups"], 4)

    def test_audit_still_requires_preds(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = _run_cli([
                "audit",
                "--matrix", str(DEMO / "m_demo.npz"),
                "--train-matrix", str(DEMO / "m_train_demo.npz"),
                "--out", str(Path(tmp) / "out"),
                "--n-boot", "10",
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--preds is required", proc.stderr)


class TestMatrixOrientation(unittest.TestCase):
    def test_npz_explicit_transpose(self):
        m = np.load(FIXTURE / "m_test_fixture.npz")
        dense = sparse.csr_matrix(
            (m["data"], m["indices"], m["indptr"]), shape=tuple(m["shape"])
        ).toarray()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            _save_csr(tmp / "t.npz", dense.T)
            got, info = io.load_matrix(
                tmp / "t.npz", orientation="cells-x-peaks", return_info=True
            )
            self.assertTrue(info["transposed"])
            np.testing.assert_array_equal(got.toarray(), dense)
            got2, info2 = io.load_matrix(tmp / "t.npz", return_info=True)
            self.assertFalse(info2["transposed"])
            np.testing.assert_array_equal(got2.toarray(), dense.T)

    def test_infer_h5ad_orientation_from_names(self):
        peaks = ["chr1:100-200", "chr2:300-400", "chr3:500-600"]
        cells = ["AAACCTGAGAAACCAT-1", "AAACCTGAGAAACCGT-1", "AAACCTGAGAAACCTA-1"]
        self.assertEqual(io.infer_h5ad_orientation(cells, peaks), "cells-x-peaks")
        self.assertEqual(io.infer_h5ad_orientation(peaks, cells), "peaks-x-cells")
        # 名字无法判断时退回 AnnData 约定
        self.assertEqual(
            io.infer_h5ad_orientation(["a", "b", "c"], ["x", "y", "z"]),
            "cells-x-peaks",
        )

    def test_h5ad_both_layouts_load_to_peaks_x_cells(self):
        try:
            import anndata as ad
            import pandas as pd
        except ImportError:  # pragma: no cover
            self.skipTest("anndata not installed")
        rng = np.random.default_rng(0)
        peaks = [f"chr1:{i * 100}-{i * 100 + 50}" for i in range(6)]
        cells = [f"AAACCTGAGAAACC{i:02d}-1" for i in range(8)]
        dense = (rng.random((6, 8)) < 0.3).astype(np.float32)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # AnnData 约定：obs x var = 细胞 x 峰
            p1 = tmp / "cells_by_peaks.h5ad"
            ad.AnnData(
                dense.T,
                obs=pd.DataFrame(index=cells),
                var=pd.DataFrame(index=peaks),
            ).write_h5ad(p1)
            # 反惯例：obs x var = 峰 x 细胞
            p2 = tmp / "peaks_by_cells.h5ad"
            ad.AnnData(
                dense,
                obs=pd.DataFrame(index=peaks),
                var=pd.DataFrame(index=cells),
            ).write_h5ad(p2)
            for path, expect_transposed in [(p1, True), (p2, False)]:
                got, info = io.load_matrix(path, return_info=True)
                self.assertEqual(info["orientation"], "cells-x-peaks"
                                 if expect_transposed else "peaks-x-cells")
                self.assertEqual(info["transposed"], expect_transposed)
                np.testing.assert_array_equal(got.toarray(), dense)
            # 显式覆盖推断
            got, info = io.load_matrix(
                p1, orientation="peaks-x-cells", return_info=True
            )
            np.testing.assert_array_equal(got.toarray(), dense.T)
            self.assertEqual(info["inferred_orientation"], "cells-x-peaks")

    def test_cli_shape_mismatch_hints_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            m = np.load(FIXTURE / "m_test_fixture.npz")
            dense = sparse.csr_matrix(
                (m["data"], m["indices"], m["indptr"]), shape=tuple(m["shape"])
            ).toarray()
            _save_csr(tmp / "t.npz", dense.T)
            proc = _run_cli([
                "depth",
                "--matrix", str(tmp / "t.npz"),
                "--preds", str(FIXTURE / "preds_fixture.npz"),
                "--out", str(tmp / "out"),
            ], cwd=tmp)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--matrix-orientation", proc.stderr)


if __name__ == "__main__":
    unittest.main()
