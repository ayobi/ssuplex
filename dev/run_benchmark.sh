#!/usr/bin/env bash
# Run SSUplex and Metaxa2 on the same reads with the same HMM profiles, then
# compare per-read origin calls.
#
# Requires on PATH: cargo (or a prebuilt binary), metaxa2, hmmsearch, python3.
#
# Usage:
#   dev/run_benchmark.sh --db <metaxa2_db> --reads <reads.fasta> --out <dir> \
#                        [--threads N] [--bin path/to/ssuplex]
#
# <metaxa2_db> is Metaxa2's database dir (contains SSU/HMMs/). On a bioconda
# install it is typically under $CONDA_PREFIX/share/metaxa2*/metaxa2_db or
# alongside the metaxa2 script — run `metaxa2 --help` or `which metaxa2` to find it.
set -euo pipefail

DB=""; READS=""; OUT=""; THREADS=4; BIN=""; TRUTH=""; MIN_DOMAINS=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)          DB="$2"; shift 2 ;;
    --reads)       READS="$2"; shift 2 ;;
    --out)         OUT="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --bin)         BIN="$2"; shift 2 ;;
    --truth)       TRUTH="$2"; shift 2 ;;
    --min-domains) MIN_DOMAINS="$2"; shift 2 ;;
    *) echo "unknown arg: $1"; exit 1 ;;
  esac
done
[[ -n "$DB" && -n "$READS" && -n "$OUT" ]] || {
  echo "usage: $0 --db <metaxa2_db> --reads <reads.fasta> --out <dir> [--threads N] [--bin BIN]"; exit 1; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in metaxa2 hmmsearch python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FAIL: '$tool' not on PATH"; exit 1; }
done

# Optional resource accounting.
TIME=""
if command -v /usr/bin/time >/dev/null 2>&1; then TIME="/usr/bin/time -v"; fi

mkdir -p "$OUT"/{hmms,metaxa2,ssuplex}

# --- binary ---------------------------------------------------------------
if [[ -z "$BIN" ]]; then
  echo ">> building SSUplex (release)"
  ( cd "$ROOT" && cargo build --release )
  BIN="$ROOT/target/release/ssuplex"
fi
[[ -x "$BIN" ]] || { echo "FAIL: binary not executable: $BIN"; exit 1; }

# --- 1. assemble HMMs from the Metaxa2 DB ---------------------------------
echo; echo "=== 1/4  assembling HMM profiles from Metaxa2 DB ==="
bash "$ROOT/dev/assemble_metaxa2_hmms.sh" "$DB" "$OUT/hmms"

# --- 2. run Metaxa2 -------------------------------------------------------
echo; echo "=== 2/4  running Metaxa2 ==="
M2_PREFIX="$OUT/metaxa2/run"
$TIME metaxa2 -i "$READS" -o "$M2_PREFIX" --cpu "$THREADS" -g ssu \
  2> "$OUT/metaxa2/metaxa2.time.log" || {
    echo "Metaxa2 run failed — see $OUT/metaxa2/metaxa2.time.log"; exit 1; }
M2_RESULTS="$M2_PREFIX.extraction.results"
[[ -f "$M2_RESULTS" ]] || { echo "FAIL: expected $M2_RESULTS not found"; ls -la "$OUT/metaxa2"; exit 1; }

# --- 3. run SSUplex ----------------------------------------------------
echo; echo "=== 3/4  running SSUplex ==="
MR_PREFIX="$OUT/ssuplex/sample"
$TIME "$BIN" -i "$READS" -o "$MR_PREFIX" --hmm-dir "$OUT/hmms" -t "$THREADS" --min-domains "$MIN_DOMAINS" \
  2> "$OUT/ssuplex/ssuplex.time.log"
MR_TSV="$MR_PREFIX.extraction.tsv"

# --- 4. compare -----------------------------------------------------------
echo; echo "=== 4/4  comparing origin calls ==="
TRUTH_ARG=()
[[ -n "$TRUTH" ]] && TRUTH_ARG=(--truth "$TRUTH")
python3 "$ROOT/dev/compare_origins.py" \
  --ssuplex "$MR_TSV" \
  --metaxa2 "$M2_RESULTS" \
  --reads "$READS" \
  "${TRUTH_ARG[@]}" \
  --out "$OUT/concordance" | tee "$OUT/concordance.txt"

echo
echo ">> artefacts in $OUT/:"
echo "   hmms/                assembled per-origin profiles"
echo "   metaxa2/run.*        Metaxa2 outputs (+ metaxa2.time.log for time/RSS)"
echo "   ssuplex/sample.*  SSUplex outputs (+ ssuplex.time.log)"
echo "   concordance.txt      the comparison printed above"
echo "   concordance.join.csv per-read join (read_id, metaxa2, ssuplex)"
