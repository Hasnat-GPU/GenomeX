#!/usr/bin/env bash
# Run real BUSCO over the demo genomes so its calls can be diffed against
# GenomeX's. Uses the same local lineage directory GenomeX uses, offline, so
# both tools score against byte-identical profiles and cutoffs.
#
#   micromamba create -y -n busco -c conda-forge -c bioconda busco=5.8.3
#   bench/run_busco.sh
#
# Two modes, because they answer different questions:
#
#   genome   -- BUSCO does its own gene prediction. End-to-end comparison: what
#               a user actually experiences choosing one tool over the other.
#   proteins -- BUSCO scores GenomeX's own predicted proteins. Isolates marker
#               classification from gene calling, so a disagreement here is
#               about cutoffs rather than about Prodigal.
#
# Running only the first would leave every difference unattributable.
#
# Env: GENOMEX_DATA, GENOMEX_LINEAGE, GENOMEX_RUN, BUSCO_ENV, BENCH_OUT, BUSCO_CPU
set -euo pipefail

DATA="${GENOMEX_DATA:-$HOME/genomex-work/data/Bacteria}"
LINEAGE="${GENOMEX_LINEAGE:-$HOME/genomex-work/db/bacteria_odb10}"
GENOMEX_RUN="${GENOMEX_RUN:-runs/demo-rhizobia}"
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
mkdir -p "$OUT/genome" "$OUT/proteins"

for g in "${GENOMES[@]}"; do
  if [ -d "$OUT/genome/$g" ]; then
    echo "== $g genome mode already done, skipping"
  else
    echo "== BUSCO $g -- genome mode"
    "$MM" run -n "$BUSCO_ENV" busco \
      -i "$DATA/$g.fna" -l "$LINEAGE" -o "$g" --out_path "$OUT/genome" \
      -m genome --offline --cpu "$CPU" -f
  fi

  faa="$GENOMEX_RUN/genomes/$g/proteins.faa"
  if [ ! -f "$faa" ]; then
    echo "   no GenomeX proteins at $faa -- run demo.sh first; skipping proteins mode"
  elif [ -d "$OUT/proteins/$g" ]; then
    echo "== $g proteins mode already done, skipping"
  else
    echo "== BUSCO $g -- proteins mode over GenomeX gene calls"
    "$MM" run -n "$BUSCO_ENV" busco \
      -i "$faa" -l "$LINEAGE" -o "$g" --out_path "$OUT/proteins" \
      -m proteins --offline --cpu "$CPU" -f
  fi
done

echo
echo "== summaries"
for mode in genome proteins; do
  echo "-- $mode"
  grep -h "C:" "$OUT/$mode"/*/short_summary.*.txt 2>/dev/null | sed 's/^[[:space:]]*//' || true
done
