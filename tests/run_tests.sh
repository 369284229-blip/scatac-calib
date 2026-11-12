#!/usr/bin/env bash
# Run the full test suite without pytest.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python -m unittest discover -s tests -v
