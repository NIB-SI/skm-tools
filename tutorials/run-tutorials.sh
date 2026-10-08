#!/usr/bin/env bash
# Run the tutorials in place, Cytoscape cells included (CI skips those), e.g. before a release.
#
# Needs a running Cytoscape (3.10 or later, with the enhancedGraphics app), and the
# tutorials extra (pip install ".[tutorials]"). To keep Cytoscape from taking the focus
# while the notebooks run, run it in Xephyr (Linux), e.g.:
#
#   Xephyr :2 -screen 1600x1000 -resizeable &
#   cd /path/to/Cytoscape && DISPLAY=:2 bash -c 'sleep infinity | ./cytoscape.sh' &
#
# Usage: ./run-tutorials.sh [notebook ...]   (default: all tutorial-*.ipynb)
set -euo pipefail
cd "$(dirname "$0")"

curl -s -m 5 http://127.0.0.1:1234/v1/version > /dev/null \
    || { echo "Cytoscape is not running (no answer on port 1234)." >&2; exit 1; }

notebooks=("$@")
[ ${#notebooks[@]} -eq 0 ] && notebooks=(tutorial-*.ipynb)
for nb in "${notebooks[@]}"; do
    echo "Running $nb"
    jupyter execute --inplace --timeout=1800 "$nb"
done
