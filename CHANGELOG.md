# Changelog

All notable changes to `scatac-calib` are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-11-12

First installable release. **Algorithm behaviour is unchanged from v0.1**;
this release only adds packaging, tests, demo data and documentation.

### Added

- `pyproject.toml` (setuptools) with the `scatac-calib` console script,
  Python ≥ 3.9, and numpy / scipy / scikit-learn dependencies.
- `tests/test_zero_model.py`: synthetic-matrix checks for the frequency-null
  model (per-peak / per-cell AUROC in 0.45–0.55 under random labels,
  Δ ≈ 0 when the model *is* the null, single-class rows/columns skipped).
- `tests/test_cli_smoke.py`: end-to-end `audit` on the bundled demo inputs.
- `tests/run_tests.sh`: `python -m unittest discover -s tests -v`.
- `examples/demo/`: 100 peaks × 80 cells synthetic inputs, a seeded
  generation script, `run_demo.sh` and `expected.md`.
- `LICENSE` (MIT), `CITATION.cff`, `CHANGELOG.md`.
- `docs/reproducibility.md`: shortest path from raw matrices to Table 2.
- `examples/real_fixture/`: fixed-index **real wheat data** subset (100 peaks ×
  80 cells from `dataset_v3` test, seed-0 predictions), `run_real_fixture.sh`,
  `expected_real.md` and `selection_manifest.json` with row indices and SHA256.
  The synthetic demo is kept and labelled as “synthetic demo” in the README.
- `docs/table3_executable.md`: every one of the 15 Table 3 rows mapped to a
  command/input/output, with `external script` paths for M1, MAPQ and
  random-peak training, plus the decision-tree “items 1–8” → 15-row mapping.
- `CODE_AVAILABILITY.md`, `DATA_AVAILABILITY.md`, `RELEASE_CHECKLIST.md`,
  `VERSION`, and a git repository tagged `v1.0.0`.  Repository URL and Zenodo
  DOI remain `{TBD}` (author decision; never fabricated).
- README: install, CLI, demo, tests, input format, Table 2/3 mapping, and
  the v1.0 regression evidence from work order 031.

### Verified

- Regression against `bootstrap_summary.json` on the wheat dataset, three
  seeds: Δ point-estimate differences are exactly 0; the largest
  confidence-interval endpoint difference is `1.51e-4`.
- 031 formal run: per-peak AUC 0.76459 ± 0.00108, per-cell AUC
  0.70465 ± 0.00170; Δ +0.06878 ± 0.00108 / +0.01684 ± 0.00170; all
  paired-bootstrap p < 0.001.

### Unchanged

- All five subcommands (`audit`, `zero-model`, `oracle`, `delta`, `depth`),
  their scoring rules, bootstrap RNG usage and output JSON schema are
  identical to v0.1.

## [0.1.0] - 2026-11-12

- Initial source-only skeleton: frequency-null model, group-level oracle,
  ΔAUROC with paired bootstrap, depth-matched per-peak AUROC, reporting
  checklist, and the `audit` pipeline.
