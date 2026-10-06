#!/usr/bin/env bash
# Accuracy and agreement with Metaxa2 on the labelled and real datasets:
#   testfasta  Metaxa2's test.fasta (50 SSU, ten per origin, plus 50 LSU negatives)
#   ref        reference set sampled from the Metaxa2 SSU database (make_truth_set.py)
#   zymo       first N ZymoBIOMICS ONT reads (SRR10391201), all bacterial
#   rice       first N rice-root ONT reads (SRR25243163), no per-read truth
# Each dataset goes through dev/run_benchmark.sh; reports are collected in
# <out>/ALL_REPORTS.txt. Re-running skips datasets that already have a report.
#
# Requires the environment described in dev/BENCHMARK.md (metaxa2, nhmmer,
# fastacmd, python3, curl).
#
# Usage:
#   dev/run_real_data.sh --metaxa2-dir <Metaxa2_2.2.3> --out <dir> \
#       [--threads 16] [--n 5000] [--bin path/to/ssuplex]
set -euo pipefail

M2DIR=""; OUT=""; THREADS=16; N=5000; BIN=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --metaxa2-dir) M2DIR="$2"; shift 2 ;;
    --out)         OUT="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --n)           N="$2"; shift 2 ;;
    --bin)         BIN="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$M2DIR" && -n "$OUT" ]] || {
  echo "usage: $0 --metaxa2-dir <Metaxa2_2.2.3> --out <dir> [--threads N] [--n READS] [--bin BIN]" >&2; exit 1; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
M2DIR="$(cd "$M2DIR" && pwd)"
DB="$M2DIR/metaxa2_db"
TF="$M2DIR/test.fasta"
[[ -d "$DB/SSU" && -s "$TF" ]] || { echo "FAIL: $M2DIR does not look like the Metaxa2 2.2.3 download" >&2; exit 1; }
if [[ -z "$BIN" ]]; then
  ( cd "$ROOT" && cargo build --release --locked --quiet )
  BIN="$ROOT/target/release/ssuplex"
fi
mkdir -p "$OUT/data"
OUT="$(cd "$OUT" && pwd)"
D="$OUT/data"

echo "=== preparing datasets ==="
n_tf="$(grep -c '^>' "$TF")"
[[ "$n_tf" -eq 100 ]] || { echo "FAIL: expected 100 records in $TF, found $n_tf" >&2; exit 1; }
awk 'BEGIN{print "read_id\ttrue_origin"; m["A"]="archaea"; m["B"]="bacteria";
           m["C"]="chloroplast"; m["E"]="eukaryota"; m["M"]="mitochondria"}
     /^>/{k++; id=substr($1,2); print id "\t" ((k<=50) ? "(negative)" : m[substr(id,length(id))])}' \
  "$TF" > "$D/testfasta.truth.tsv"

[[ -s "$D/ref.reads.fasta" ]] || python3 "$ROOT/dev/make_truth_set.py" --db "$DB" --out "$D/ref" \
  --per-origin 100 --negatives 10 --seed 42

[[ -s "$D/zymo.fasta" ]] || "$ROOT/dev/fetch_ena_reads.sh" SRR10391201 "$N" "$D/zymo.fasta"
awk 'BEGIN{print "read_id\ttrue_origin"} /^>/{print substr($1,2) "\tbacteria"}' \
  "$D/zymo.fasta" > "$D/zymo.truth.tsv"

[[ -s "$D/rice.fasta" ]] || "$ROOT/dev/fetch_ena_reads.sh" SRR25243163 "$N" "$D/rice.fasta"

run() {  # run <name> <reads> [truth]
  local name="$1" reads="$2" truth="${3:-}"
  if [[ -s "$OUT/$name/report.txt" ]]; then
    echo "=== $name: report exists, skipping ==="; return 0
  fi
  echo; echo "=== $name ==="
  local targs=()
  [[ -n "$truth" ]] && targs=(--truth "$truth")
  "$ROOT/dev/run_benchmark.sh" --db "$DB" --reads "$reads" --out "$OUT/$name" \
    --threads "$THREADS" --bin "$BIN" "${targs[@]+"${targs[@]}"}"
}

run testfasta "$TF" "$D/testfasta.truth.tsv"
run ref "$D/ref.reads.fasta" "$D/ref.truth.tsv"
run zymo "$D/zymo.fasta" "$D/zymo.truth.tsv"
run rice "$D/rice.fasta"

cat "$OUT/testfasta/versions.txt" > "$OUT/ALL_REPORTS.txt"
for name in testfasta ref zymo rice; do
  echo >> "$OUT/ALL_REPORTS.txt"
  cat "$OUT/$name/report.txt" >> "$OUT/ALL_REPORTS.txt"
done
echo; echo "Done. Combined report: $OUT/ALL_REPORTS.txt"
