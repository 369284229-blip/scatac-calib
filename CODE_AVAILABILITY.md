# Code Availability — scatac-calib

> Status: release candidate prepared by worker 050 (2026-11-12).  Items marked
> `{TBD}` are decisions only the authors can make; they have **not** been
> invented or filled with a fake URL/DOI.

## Software identity

| Field | Value |
|---|---|
| Name | `scatac-calib` (`scatac_calib` import name) |
| Version | **1.0.0** (`VERSION`, `scatac_calib.__version__`, `pyproject.toml`) |
| Release date | 2026-11-12 |
| Git tag | **`v1.0.0`** |
| Commit | tag `v1.0.0` target; the exact hash is in the shipped **`COMMIT`** file (`git rev-parse HEAD` at packaging time) and can be recovered with `git rev-parse v1.0.0` |
| License | **MIT** (`LICENSE`; SPDX `MIT`) |
| Copyright holder | `Wheat scATAC Benchmark Team` (placeholder — real holder `{TBD}`) |
| Public repository URL | **`{TBD}`** (user must decide public/private and host) |
| Archived DOI (Zenodo) | **`{TBD}`** |
| Author list / affiliation | `{TBD}` (`CITATION.cff` currently uses the team placeholder) |
| Contact | `{TBD}` |

## What the release contains

- Five subcommands: `audit`, `zero-model`, `oracle`, `delta`, `depth`.
- Packaged data: a **synthetic demo** (`examples/demo/`) and a **real-data
  fixture** (`examples/real_fixture/`, 100 peaks × 80 cells from the wheat
  `dataset_v3` test split; provenance in
  `examples/real_fixture/selection_manifest.json`).
- Tests that do not require pytest (`tests/`), the decision tree
  (`docs/decision_tree.md`), the Table 3 executable mapping
  (`docs/table3_executable.md`) and the shortest reproduction path
  (`docs/reproducibility.md`).
- No model weights, no raw sequencing data, no full cell × peak matrices.

## Installation

From the wheel (no network needed if dependencies are present):

```bash
python -m pip install dist/scatac_calib-1.0.0-py3-none-any.whl
scatac-calib --help
```

From the sdist:

```bash
python -m pip install dist/scatac_calib-1.0.0.tar.gz
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
# 1) unit tests + CLI smoke test (9 tests)
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

> The `scatac-calib` toolkit (v1.0.0, tag `v1.0.0`, commit hash in the
> accompanying `COMMIT` file) is released under the MIT licence at
> `{REPOSITORY_URL}` and archived at Zenodo under DOI `{ZENODO_DOI}`.  Installation, the synthetic demo, the real-data
> fixture, tests and the mapping from the reporting checklist to CLI commands
> are included in the archive.

Do **not** paste this into the manuscript until the URL and DOI exist.
`写作/Methods_draft.md §6.10` currently says “v0.1” and `{TBD}` for the commit;
that is a manuscript-side fix and is out of scope for this work order.

## Version note

The package is **1.0.0** everywhere (pyproject, `__init__`, `VERSION`,
`CITATION.cff`, wheel/sdist metadata).  Older text that says “v0.1” refers to
the pre-packaging skeleton, not this release; `CHANGELOG.md` records both.
