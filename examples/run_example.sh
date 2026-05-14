#!/usr/bin/env bash
# Example: run the pipeline on every .rnx file in a directory.

set -euo pipefail

INPUT_DIR="${1:-./examples/data}"
OUT_DIR="${2:-./out}"

mkdir -p "$OUT_DIR"
for f in "$INPUT_DIR"/*.rnx "$INPUT_DIR"/*.obs; do
  [ -e "$f" ] || continue
  stem=$(basename "$f" | sed 's/\.[^.]*$//')
  echo "── $stem ──"
  python -m gnss_nav.cli "$f" \
    --out-dir "$OUT_DIR/$stem" \
    --elevation-mask 10 \
    --systems G,E \
    --rate 1
done
