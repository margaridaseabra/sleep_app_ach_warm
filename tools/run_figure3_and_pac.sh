#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
python tools/build_figure3_and_pac.py --repo "$PWD"
