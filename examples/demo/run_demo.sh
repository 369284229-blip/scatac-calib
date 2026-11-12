#!/usr/bin/env bash
# One-command demo: audit the bundled 100-peak x 80-cell synthetic dataset.
#
# Usage:
#   ./run_demo.sh [OUTDIR]      # default: examples/demo/demo_out
#
# Requires scatac-calib to be installed (`pip install .`) or run from the
# repository root.  Falls back to `python -m scatac_calib` when the console
# script is not on PATH.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$HERE/demo_out}"

if command -v scatac-calib >/dev/null 2>&1; then
  CMD=(scatac-calib)
else
  CMD=(env "PYTHONPATH=$HERE/../..${PYTHONPATH:+:$PYTHONPATH}" python -m scatac_calib)
fi

cd "$HERE"
"${CMD[@]}" audit \
  --matrix m_demo.npz \
  --preds preds_demo.npz \
  --train-matrix m_train_demo.npz \
  --groups groups_demo.txt \
  --out "$OUT" \
  --n-boot 200 \
  --depth-bins 5

echo "wrote $OUT/calib_report.json and $OUT/reporting_checklist.md"
