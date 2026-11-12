"""End-to-end smoke test: `audit` on the bundled demo inputs."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "examples" / "demo"


class TestAuditDemo(unittest.TestCase):
    def test_audit_writes_report_and_checklist(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "demo_out"
            cmd = [
                sys.executable, "-m", "scatac_calib", "audit",
                "--matrix", str(DEMO / "m_demo.npz"),
                "--preds", str(DEMO / "preds_demo.npz"),
                "--train-matrix", str(DEMO / "m_train_demo.npz"),
                "--groups", str(DEMO / "groups_demo.txt"),
                "--out", str(out),
                "--n-boot", "50",
                "--depth-bins", "5",
            ]
            env = dict(os.environ)
            env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
            proc = subprocess.run(
                cmd, cwd=tmp, env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)

            report_path = out / "calib_report.json"
            checklist_path = out / "reporting_checklist.md"
            self.assertTrue(report_path.exists())
            self.assertTrue(checklist_path.exists())

            result = json.loads(report_path.read_text())
            self.assertEqual(result["n_peaks"], 100)
            self.assertEqual(result["n_cells"], 80)
            for key in ("metrics", "zero_model", "delta_vs_zero", "ci95",
                        "oracle", "depth_matched", "ap"):
                self.assertIn(key, result)

            self.assertGreater(result["metrics"]["per_peak_auc"], 0.5)
            self.assertGreater(result["metrics"]["per_cell_auc"], 0.5)
            self.assertGreater(
                result["delta_vs_zero"]["peak"]["delta"], 0.0,
                msg="demo predictions should beat the frequency null",
            )
            self.assertEqual(
                result["zero_model"]["cell_freq_source"], "train",
            )
            self.assertEqual(result["oracle"]["source"], "train")
            self.assertEqual(result["oracle"]["n_groups"], 4)
            self.assertEqual(result["depth_matched"]["n_bins"], 5)


if __name__ == "__main__":
    unittest.main()
