"""End-to-end test: `audit` on the bundled real-data fixture.

The fixture is a fixed 100-peak x 80-cell subset of the wheat dataset_v3 test
split (see examples/real_fixture/selection_manifest.json).  Point estimates do
not depend on --n-boot, so a small bootstrap count keeps the test fast; the
reference values are in examples/real_fixture/expected_real.md.
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "examples" / "real_fixture"


class TestAuditRealFixture(unittest.TestCase):
    def test_audit_matches_expected_point_estimates(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "fixture_out"
            cmd = [
                sys.executable, "-m", "scatac_calib", "audit",
                "--matrix", str(FIXTURE / "m_test_fixture.npz"),
                "--preds", str(FIXTURE / "preds_fixture.npz"),
                "--train-matrix", str(FIXTURE / "m_train_fixture.npz"),
                "--groups", str(FIXTURE / "groups_fixture.txt"),
                "--out", str(out),
                "--n-boot", "50",
                "--depth-bins", "5",
                "--seed", "20261112",
            ]
            env = dict(os.environ)
            env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
            proc = subprocess.run(
                cmd, cwd=tmp, env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)

            r = json.loads((out / "calib_report.json").read_text())
            self.assertEqual(r["n_peaks"], 100)
            self.assertEqual(r["n_cells"], 80)
            self.assertEqual(r["zero_model"]["cell_freq_source"], "train")
            self.assertEqual(r["oracle"]["source"], "test(label-informed)")
            self.assertEqual(r["oracle"]["n_groups"], 14)
            self.assertEqual(r["metrics"]["n_valid_peaks"], 73)
            self.assertEqual(r["metrics"]["n_valid_cells"], 64)
            self.assertEqual(r["depth_matched"]["n_bins"], 5)

            expected = {
                ("metrics", "per_peak_auc"): 0.747824,
                ("metrics", "per_cell_auc"): 0.694087,
                ("ap",): 0.077476,
                ("zero_model", "per_peak_auc"): 0.683415,
                ("zero_model", "per_cell_auc"): 0.808941,
                ("oracle", "per_peak_auc"): 0.934845,
                ("oracle", "per_cell_auc"): 0.955008,
                ("depth_matched", "per_peak_auc"): 0.623545,
                ("delta_vs_zero", "peak", "delta"): 0.064409,
                ("delta_vs_zero", "cell", "delta"): -0.114853,
            }
            for path, want in expected.items():
                got = r
                for key in path:
                    got = got[key]
                self.assertAlmostEqual(
                    got, want, delta=1e-6,
                    msg=f"{'.'.join(path)} = {got}, expected {want}",
                )


if __name__ == "__main__":
    unittest.main()
