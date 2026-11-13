# Code Availability — scatac-calib

> Status: release candidate prepared by worker 050 (2026-11-12).  Items marked
> `{TBD}` are decisions only the authors can make; they have **not** been
> invented or filled with a fake URL/DOI.

## Software identity

| Field | Value |
|---|---|
| Name | `scatac-calib` (`scatac_calib` import name) |
| Version | **1.1.0** (`VERSION`, `scatac_calib.__version__`, `pyproject.toml`) |
| Release date | 2026-11-12 |
| Git tag | **`v1.0.0`** for the stage-1 release; `v1.1.0` tag not created yet (worker 068 does not push/tag) |
| Commit | v1.0.0 tag target `a6d457f` (the shipped **`COMMIT`** file); the v1.1.0 working tree is **not committed** by worker 068, so the commander must commit and regenerate `COMMIT` before tagging `v1.1.0` |
| License | **MIT** (`LICENSE`; SPDX `MIT`) |
| Copyright holder | `Wheat scATAC Benchmark Team` (placeholder — real holder `{TBD}`) |
| Public repository URL | **`https://github.com/369284229-blip/scatac-calib`** (public; `v1.0.0` tagged, `v1.1.0` pending) |
| Archived DOI (Zenodo) | **`{TBD}`** |
| Author list / affiliation | `{TBD}` (`CITATION.cff` currently uses the team placeholder) |
| Contact | `{TBD}` |

## What the release contains

- Five subcommands: `audit`, `zero-model`, `oracle`, `delta`, `depth`.
- Packaged data: a **synthetic demo** (`examples/demo/`) and a **real-data
  fixture** (`examples/real_fixture/`, 100 peaks × 80 cells from the wheat
  `dataset_v3` test split; provenance in
  `examples/real_fixture/selection_manifest.json`).
- `--null inductive` (5-fold cell cross-fitting) and input hard validation
  (`--train-matrix`, `.h5ad` orientation, optional `--preds` for zero/oracle).
- `examples/full_reproduction/reproduce_inductive.py`, a read-only script that
  reproduces the full-scale null 0.6852978 / Δcell +0.0193566 from project
  data.  The 430 MB of full matrices/predictions are not packaged.
- Tests that do not require pytest (`tests/`), the decision tree
  (`docs/decision_tree.md`), the Table 3 executable mapping
  (`docs/table3_executable.md`) and the shortest reproduction path
  (`docs/reproducibility.md`).
- No model weights, no raw sequencing data, no full cell × peak matrices.

## Installation

From the wheel (no network needed if dependencies are present):

```bash
python -m pip install dist/scatac_calib-1.1.0-py3-none-any.whl
scatac-calib --help
```

From the sdist:

```bash
python -m pip install dist/scatac_calib-1.1.0.tar.gz
```

From a source checkout (development):

```bash
python -m pip install -e .
scatac-calib --help
```

Requirements: Python ≥ 3.9, `numpy`, `scipy`, `scikit-learn`.  Reading `.h5ad`
additionally needs `scanpy` (`pip install "scatac-calib[h5ad]"`).

## Verification (what a third party can run)

```bash
# 1) unit tests + CLI tests
python -m unittest discover -s tests -v

# 2) synthetic demo (no data access needed)
cd examples/demo && ./run_demo.sh

# 3) real-data fixture (shipped subset; real wheat data)
cd examples/real_fixture && ./run_real_fixture.sh
```

Expected numbers and tolerances: `examples/demo/expected.md` and
`examples/real_fixture/expected_real.md` (point estimates 1e-6, bootstrap CIs
1e-3).  The v1.0 regression against the project results is documented in
`README.md` §“已验证”.

## Suggested Code Availability statement (fill the two placeholders)

> The `scatac-calib` toolkit (v1.1.0, tag `v1.1.0`, commit hash in the
> accompanying `COMMIT` file) is released under the MIT licence at
> `https://github.com/369284229-blip/scatac-calib` and archived at Zenodo under DOI `{ZENODO_DOI}`.  Installation, the synthetic demo, the real-data
> fixture, tests and the mapping from the reporting checklist to CLI commands
> are included in the archive.

Do **not** paste this into the manuscript until the DOI exists.
`写作/Methods_draft.md §6.10` currently says “v0.1” and `{TBD}` for the commit;
that is a manuscript-side fix and is out of scope for this work order.

## Version note

The package is **1.1.0** everywhere (pyproject, `__init__`, `VERSION`,
`CITATION.cff`, wheel/sdist metadata).  v1.0.0 remains the tag/DOI line for the
stage-1 release; v1.1.0 is the stage-2 tool update (inductive null + input
validation).  Older text that says “v0.1” refers to the pre-packaging skeleton;
`CHANGELOG.md` records all three.
