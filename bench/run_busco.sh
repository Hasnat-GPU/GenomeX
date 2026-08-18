#!/usr/bin/env bash
# Run real BUSCO over the demo genomes so its calls can be diffed against
# GenomeX's. Uses the same local lineage directory GenomeX uses, offline, so
# both tools score against byte-identical profiles and cutoffs.
#
#   micromamba create -y -n busco -c conda-forge -c bioconda busco
#   bench/run_busco.sh
#
# Env: GENOMEX_DATA, GENOMEX_LINEAGE, BUSCO_ENV, BENCH_OUT, BUSCO_CPU
set -euo pipefail

DATA="${GENOMEX_DATA:-$HOME/genomex-work/data/Bacteria}"
LINEAGE="${GENOMEX_LINEAGE:-$HOME/genomex-work/db/bacteria_odb10}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
OUT="${BENCH_OUT:-bench/busco_out}"
CPU="${BUSCO_CPU:-$(nproc 2>/dev/null || echo 4)}"
MM="${MICROMAMBA:-$HOME/.local/bin/micromamba}"

GENOMES=(
  GCA_000300095.1
  GCA_040948545.1
  GCF_000020045.1
  GCA_000069785.1
)

test -d "$LINEAGE/hmms" || { echo "no lineage at $LINEAGE" >&2; exit 1; }
mkdir -p "$OUT"

for g in "${GENOMES[@]}"; do
  if [ -d "$OUT/$g" ]; then
    echo "== $g already done, skipping"
    continue
  fi
  echo "== BUSCO $g"
  "$MM" run -n "$BUSCO_ENV" busco \
    -i "$DATA/$g.fna" \
    -l "$LINEAGE" \
    -o "$g" \
    --out_path "$OUT" \
    -m genome \
    --offline \
    --cpu "$CPU" \
    -f
done

echo
echo "== summaries"
grep -h "C:" "$OUT"/*/short_summary.*.txt 2>/dev/null | sed 's/^\s*//' || true
