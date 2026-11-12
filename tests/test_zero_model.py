"""Unit tests for the frequency-null (zero) model.

The checks use fully synthetic random binary matrices, so they assert
properties of the implementation rather than properties of any dataset.

Run with either of:
    python -m unittest discover -s tests -v
    ./tests/run_tests.sh
"""

import unittest

import numpy as np

from scatac_calib import baselines, bootstrap, metrics


N_PEAKS = 200
# 400 rather than 100 cells: peak_freq is estimated from the *same* test
# matrix, which gives the per-cell null a finite-sample upward bias of about
# +0.056 at 100 cells (1/sqrt(n_cells * p(1-p)) scaling).  At 400 cells the
# bias is < 0.03, so a 0.45-0.55 window is a meaningful assertion.  The
# 200-peak dimension keeps the per-peak bias negligible (~0.04 in the
# no-train fallback); the behaviour itself does not depend on matrix size.
N_CELLS = 400


def random_binary(seed, shape=(N_PEAKS, N_CELLS)):
    rng = np.random.default_rng(seed)
    return rng.integers(0, 2, size=shape).astype(np.float32)


class TestNullModelAuc(unittest.TestCase):
    """Random labels / random scores -> AUROC near 0.5 on both axes."""

    def setUp(self):
        self.labels = random_binary(11)
        self.train = random_binary(12)  # independent draw, same shape

    def test_zero_model_auc_is_near_half(self):
        peak_freq, cell_freq, source = baselines.zero_scores(self.train, self.labels)
        self.assertEqual(source, "train")
        z_peak, z_cell = baselines.zero_auc_arrays(self.labels, peak_freq, cell_freq)

        self.assertEqual(z_peak.shape, (N_PEAKS,))
        self.assertEqual(z_cell.shape, (N_CELLS,))
        self.assertTrue(np.all(np.isfinite(z_peak)))
        self.assertTrue(np.all(np.isfinite(z_cell)))

        # Every random binary row/column has both classes with overwhelming
        # probability, so the means are over all 200 peaks / 100 cells.
        self.assertGreaterEqual(float(np.mean(z_peak)), 0.45)
        self.assertLessEqual(float(np.mean(z_peak)), 0.55)
        self.assertGreaterEqual(float(np.mean(z_cell)), 0.45)
        self.assertLessEqual(float(np.mean(z_cell)), 0.55)

    def test_zero_model_without_train_uses_test_marginals(self):
        peak_freq, cell_freq, source = baselines.zero_scores(None, self.labels)
        self.assertEqual(source, "test")
        np.testing.assert_allclose(peak_freq, self.labels.mean(axis=1), rtol=0, atol=1e-7)
        np.testing.assert_allclose(cell_freq, self.labels.mean(axis=0), rtol=0, atol=1e-7)

        z_peak, z_cell = baselines.zero_auc_arrays(self.labels, peak_freq, cell_freq)
        self.assertAlmostEqual(float(np.mean(z_peak)), 0.5, delta=0.05)
        self.assertAlmostEqual(float(np.mean(z_cell)), 0.5, delta=0.05)


class TestDeltaWhenModelIsNull(unittest.TestCase):
    """If the model prediction *is* the null score, delta must be zero."""

    def setUp(self):
        self.labels = random_binary(21)
        self.train = random_binary(22)
        self.peak_freq, self.cell_freq, _ = baselines.zero_scores(self.train, self.labels)
        # Row-wise ranking follows cell_freq, column-wise ranking follows
        # peak_freq, so this prediction reproduces both null axes exactly.
        self.preds = (
            self.peak_freq[:, None] + self.cell_freq[None, :]
        ).astype(np.float32)
        self.rng = np.random.default_rng(23)

    def test_delta_is_zero_on_both_axes(self):
        model_peak = metrics.per_peak_auc(self.labels, self.preds)
        model_cell = metrics.per_cell_auc(self.labels, self.preds)
        z_peak, z_cell = baselines.zero_auc_arrays(
            self.labels, self.peak_freq, self.cell_freq
        )

        # The null scores are reproduced exactly up to a monotone constant.
        np.testing.assert_allclose(model_peak, z_peak, rtol=0, atol=1e-9)
        np.testing.assert_allclose(model_cell, z_cell, rtol=0, atol=1e-9)

        for model_auc, null_auc in ((model_peak, z_peak), (model_cell, z_cell)):
            summary = bootstrap.delta_summary(model_auc, null_auc, self.rng, n_boot=300)
            self.assertEqual(summary["delta"], 0.0)
            self.assertEqual(summary["paired"]["mean"], 0.0)
            self.assertEqual(summary["paired"]["ci95"], [0.0, 0.0])
            self.assertEqual(summary["paired"]["p"], 1.0)
            self.assertEqual(summary["paired"]["n"], len(model_auc))


class TestSingleClassSkipped(unittest.TestCase):
    """All-zero / all-one peaks and cells must produce NaN, not an error."""

    def setUp(self):
        # A dead peak and a dead cell are compatible (their intersection is 0);
        # a saturated peak and a saturated cell are compatible (intersection
        # 1).  Mixing dead and saturated would force a 0/1 contradiction at
        # the intersection, so the two cases use separate matrices.
        dead = random_binary(31)
        dead[0, :] = 0.0     # all-zero peak
        dead[:, 2] = 0.0     # all-zero cell
        self.dead = dead

        sat = random_binary(31)
        sat[5, :] = 1.0      # all-one peak
        sat[:, 7] = 1.0      # all-one cell
        self.sat = sat

        rng = np.random.default_rng(32)
        self.preds = rng.random((N_PEAKS, N_CELLS)).astype(np.float32)
        self.train = random_binary(33)

    def test_per_peak_and_per_cell_skip_all_zero(self):
        peak_auc = metrics.per_peak_auc(self.dead, self.preds)
        cell_auc = metrics.per_cell_auc(self.dead, self.preds)

        self.assertTrue(np.isnan(peak_auc[0]))
        self.assertFalse(np.isnan(peak_auc[1]))
        self.assertTrue(np.isnan(cell_auc[2]))
        self.assertFalse(np.isnan(cell_auc[1]))

        summary = metrics.mean_auc(peak_auc, cell_auc)
        self.assertEqual(summary["n_valid_peaks"], N_PEAKS - 1)
        self.assertEqual(summary["n_valid_cells"], N_CELLS - 1)
        self.assertAlmostEqual(
            summary["per_peak_auc"], float(np.nanmean(peak_auc)), places=12
        )

    def test_per_peak_and_per_cell_skip_all_one(self):
        peak_auc = metrics.per_peak_auc(self.sat, self.preds)
        cell_auc = metrics.per_cell_auc(self.sat, self.preds)

        self.assertTrue(np.isnan(peak_auc[5]))
        self.assertFalse(np.isnan(peak_auc[1]))
        self.assertTrue(np.isnan(cell_auc[7]))
        self.assertFalse(np.isnan(cell_auc[1]))

        summary = metrics.mean_auc(peak_auc, cell_auc)
        self.assertEqual(summary["n_valid_peaks"], N_PEAKS - 1)
        self.assertEqual(summary["n_valid_cells"], N_CELLS - 1)

    def test_zero_model_skips_single_class(self):
        for matrix, bad_peak, bad_cell in (
            (self.dead, 0, 2),
            (self.sat, 5, 7),
        ):
            peak_freq, cell_freq, _ = baselines.zero_scores(self.train, matrix)
            z_peak, z_cell = baselines.zero_auc_arrays(matrix, peak_freq, cell_freq)

            self.assertTrue(np.isnan(z_peak[bad_peak]))
            self.assertTrue(np.isnan(z_cell[bad_cell]))
            self.assertEqual(int(np.sum(np.isfinite(z_peak))), N_PEAKS - 1)
            self.assertEqual(int(np.sum(np.isfinite(z_cell))), N_CELLS - 1)

    def test_delta_ignores_nan_pairs(self):
        peak_auc = metrics.per_peak_auc(self.dead, self.preds)
        peak_freq, cell_freq, _ = baselines.zero_scores(self.train, self.dead)
        z_peak, _ = baselines.zero_auc_arrays(self.dead, peak_freq, cell_freq)
        rng = np.random.default_rng(34)
        summary = bootstrap.delta_summary(peak_auc, z_peak, rng, n_boot=100)
        self.assertEqual(summary["paired"]["n"], N_PEAKS - 1)
        self.assertTrue(np.isfinite(summary["delta"]))

    def test_safe_auc_returns_nan_for_single_class(self):
        self.assertTrue(np.isnan(metrics.safe_auc(np.zeros(10), np.arange(10))))
        self.assertTrue(np.isnan(metrics.safe_auc(np.ones(10), np.arange(10))))
        self.assertAlmostEqual(metrics.safe_auc(np.array([0, 0, 1, 1]),
                                                np.array([0., 1., 2., 3.])), 1.0)


if __name__ == "__main__":
    unittest.main()
