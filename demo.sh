#!/usr/bin/env bash
# Reproduce the demo run: two same-species isolates, one congener, one out-group.
# Point GENOMEX_DATA at wherever the Bacteria/ assemblies live.
set -euo pipefail
DATA="${GENOMEX_DATA:-$HOME/genomex-work/data/Bacteria}"

if [ ! -d "$DATA" ]; then
  echo "dataset not found at $DATA -- set GENOMEX_DATA to the Bacteria/ directory" >&2
  exit 1
fi

python -m genomex run \
  "$DATA/GCA_000300095.1.fna" \
  "$DATA/GCA_040948545.1.fna" \
  "$DATA/GCF_000020045.1.fna" \
  "$DATA/GCA_000069785.1.fna" \
  --outdir runs/demo-rhizobia --all-pairs
