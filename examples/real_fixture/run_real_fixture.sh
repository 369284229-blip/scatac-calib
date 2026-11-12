#!/usr/bin/env bash
# One-command audit of the bundled real-data fixture.
#
# The fixture is a fixed-index subset of the wheat dataset_v3 test split:
# 100 peaks (chr3) x 80 cells (14 Leiden clusters).  It is REAL data, not
# synthetic; provenance and indices are in selection_manifest.json.
#
# Usage:
#   ./run_real_fixture.sh [OUTDIR]      # default: examples/real_fixture/fixture_out
#
# Requires scatac-calib to be installed (`pip install .`) or run from the
# repository root.  Falls back to `python -m scatac_calib` when the console
# script is not on PATH.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$HERE/fixture_out}"

if command -v scatac-calib >/dev/null 2>&1; then
  CMD=(scatac-calib)
else
  CMD=(env "PYTHONPATH=$HERE/../..${PYTHONPATH:+:$PYTHONPATH}" python -m scatac_calib)
fi

cd "$HERE"
"${CMD[@]}" audit \
  --matrix m_test_fixture.npz \
  --preds preds_fixture.npz \
  --train-matrix m_train_fixture.npz \
  --groups groups_fixture.txt \
  --out "$OUT" \
  --n-boot 1000 \
  --depth-bins 5 \
  --seed 20261112

echo "wrote $OUT/calib_report.json and $OUT/reporting_checklist.md"
