#!/usr/bin/env bash
# Accuracy harness: run Metaxa2 and SSUplex on the same reads with the same
# HMM profiles, then score both with dev/score_origins.py.
#
# Metaxa2 runs its full default pipeline by default, which yields both its
# HMM-stage origin (.extraction.results) and its BLAST-assisted final origin
# (per-origin FASTAs). Use --m2-mode extract to run extraction only (-x T).
# SSUplex runs once per ranking statistic (mean, sum, count); the mean run also
# writes the --debug-scores tables used for the disagreement breakdown.
#
# Requires on PATH: metaxa2 (Bioconda package "metaxa"), nhmmer, python3.
#   conda create -n metaxa -c conda-forge -c bioconda metaxa
#
# Usage:
#   dev/run_benchmark.sh --db <metaxa2_db> --reads <reads.fasta> --out <dir> \
#       [--truth <truth.tsv>] [--threads N] [--bin path/to/ssuplex] \
#       [--m2-mode full|extract]
#
# <metaxa2_db> is the database directory of a Metaxa2 2.2.3 download or
# install (it contains SSU/HMMs/). The truth TSV format is described in
# dev/score_origins.py.
set -euo pipefail

DB=""; READS=""; OUT=""; THREADS=4; BIN=""; TRUTH=""; M2MODE="full"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)      DB="$2"; shift 2 ;;
    --reads)   READS="$2"; shift 2 ;;
    --out)     OUT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --bin)     BIN="$2"; shift 2 ;;
    --truth)   TRUTH="$2"; shift 2 ;;
    --m2-mode) M2MODE="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$DB" && -n "$READS" && -n "$OUT" ]] || {
  echo "usage: $0 --db <metaxa2_db> --reads <reads.fasta> --out <dir> [--truth T] [--threads N] [--bin BIN] [--m2-mode full|extract]" >&2
  exit 1; }
[[ "$M2MODE" == full || "$M2MODE" == extract ]] || { echo "--m2-mode must be full or extract" >&2; exit 1; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in metaxa2 nhmmer python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FAIL: '$tool' not on PATH" >&2; exit 1; }
done
if [[ -z "$BIN" ]]; then
  ( cd "$ROOT" && cargo build --release --locked --quiet )
  BIN="$ROOT/target/release/ssuplex"
fi
[[ -x "$BIN" ]] || { echo "FAIL: SSUplex binary not executable: $BIN" >&2; exit 1; }

mkdir -p "$OUT"/{hmms,metaxa2,ssuplex}
{
  echo "date:     $(date)"
  echo "machine:  $(uname -srm)"
  echo "ssuplex:  $("$BIN" --version)  ($BIN)"
  echo "nhmmer:   $(command -v nhmmer)  $(nhmmer -h | sed -n 2p)"
  echo "metaxa2:  $(command -v metaxa2)  $(metaxa2 --help 2>&1 | grep -i -m 1 'version' || true)"
  echo "hmmsearch used by Metaxa2: $(command -v hmmsearch)  $(hmmsearch -h | sed -n 2p)"
  echo "Metaxa2 mode: $M2MODE   threads: $THREADS   reads: $READS"
} > "$OUT/versions.txt"
cat "$OUT/versions.txt"

echo; echo "=== 1/4 assembling HMM profiles from the Metaxa2 database ==="
bash "$ROOT/dev/assemble_metaxa2_hmms.sh" "$DB" "$OUT/hmms"

echo; echo "=== 2/4 running Metaxa2 ($M2MODE) ==="
M2_ARGS=(-i "$READS" -o "$OUT/metaxa2/run" --cpu "$THREADS" -g ssu)
[[ "$M2MODE" == extract ]] && M2_ARGS+=(-x T)
metaxa2 "${M2_ARGS[@]}" > "$OUT/metaxa2/metaxa2.log" 2>&1 \
  || { echo "FAIL: Metaxa2 failed; see $OUT/metaxa2/metaxa2.log" >&2; exit 1; }

echo; echo "=== 3/4 running SSUplex (mean, sum, count) ==="
for rank in mean sum count; do
  extra=()
  [[ "$rank" == mean ]] && extra=(--debug-scores)
  "$BIN" -i "$READS" -o "$OUT/ssuplex/$rank" --hmm-dir "$OUT/hmms" -t "$THREADS" \
    --rank "$rank" --force "${extra[@]+"${extra[@]}"}" > "$OUT/ssuplex/$rank.log" 2>&1 \
    || { echo "FAIL: SSUplex ($rank) failed; see $OUT/ssuplex/$rank.log" >&2; exit 1; }
done

echo; echo "=== 4/4 scoring ==="
TRUTH_ARG=()
[[ -n "$TRUTH" ]] && TRUTH_ARG=(--truth "$TRUTH")
python3 "$ROOT/dev/score_origins.py" --name "$(basename "$READS")" --reads "$READS" \
  --metaxa2-prefix "$OUT/metaxa2/run" --ssuplex-dir "$OUT/ssuplex" --out "$OUT" \
  "${TRUTH_ARG[@]+"${TRUTH_ARG[@]}"}" | tee "$OUT/report.txt"
