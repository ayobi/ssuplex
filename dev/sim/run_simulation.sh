#!/usr/bin/env bash
# Independent simulated benchmark, end to end:
#   1. download SILVA 138.2 SSU and the RefSeq organelle genomes (fetch_sources.sh)
#   2. build labelled sources independent of the Metaxa2 database (build_sources.py)
#   3. simulate full-length reads at three accuracy levels       (simulate_reads.py)
#   4. run Metaxa2 and SSUplex on each set and score them        (../run_benchmark.sh)
#   5. break accuracy down by novelty band and subgroup          (score_by_band.py)
# Each step is skipped when its output exists, so the script can be re-run.
#
# Requires the benchmark environment (dev/BENCHMARK.md: metaxa2, nhmmer, fastacmd)
# active, plus a second environment for the simulation tools:
#   conda create -y -n ssuplex-sim -c conda-forge -c bioconda python=3.12 badread vsearch
#
# Usage:
#   dev/sim/run_simulation.sh --metaxa2-dir <Metaxa2_2.2.3> --out <dir> [--threads 16] \
#       [--per-origin 2000] [--negatives 2000] [--reads-per-source 2] \
#       [--profiles ont_r9,ont_r10,hifi_like] [--sim-env ssuplex-sim]
set -euo pipefail

M2DIR=""; OUT=""; THREADS=16; PER=2000; NEG=2000; RPS=2
PROFILES="ont_r9,ont_r10,hifi_like"; SIMENV="ssuplex-sim"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --metaxa2-dir)      M2DIR="$2"; shift 2 ;;
    --out)              OUT="$2"; shift 2 ;;
    --threads)          THREADS="$2"; shift 2 ;;
    --per-origin)       PER="$2"; shift 2 ;;
    --negatives)        NEG="$2"; shift 2 ;;
    --reads-per-source) RPS="$2"; shift 2 ;;
    --profiles)         PROFILES="$2"; shift 2 ;;
    --sim-env)          SIMENV="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$M2DIR" && -n "$OUT" ]] || { echo "usage: $0 --metaxa2-dir <Metaxa2_2.2.3> --out <dir> [options]" >&2; exit 1; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DB="$(cd "$M2DIR" && pwd)/metaxa2_db"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
for tool in metaxa2 nhmmer fastacmd python3 curl; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FAIL: '$tool' not on PATH (activate the benchmark environment)" >&2; exit 1; }
done

# 'conda' is usually a shell function, which scripts do not inherit; use the
# executable that conda records in CONDA_EXE when an environment is activated.
CONDA_BIN="${CONDA_EXE:-$(command -v conda || true)}"
[[ -n "$CONDA_BIN" && -x "$CONDA_BIN" ]] || {
  echo "FAIL: cannot find the conda executable (CONDA_EXE is not set and conda is not on PATH)" >&2; exit 1; }
SIM=("$CONDA_BIN" run --no-capture-output -n "$SIMENV")
"$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$SIMENV" || {
  echo "FAIL: conda environment '$SIMENV' not found. Create it with:" >&2
  echo "  conda create -y -n $SIMENV -c conda-forge -c bioconda python=3.12 badread vsearch" >&2
  exit 1; }
if ! check="$("${SIM[@]}" python3 -c "import badread, edlib" 2>&1)"; then
  echo "FAIL: badread or edlib cannot be imported in '$SIMENV':" >&2; echo "$check" >&2; exit 1
fi
"${SIM[@]}" vsearch --version >/dev/null 2>&1 || { echo "FAIL: vsearch does not run in '$SIMENV'" >&2; exit 1; }
( cd "$ROOT" && cargo build --release --locked --quiet )
BIN="$ROOT/target/release/ssuplex"

echo "=== 1/5 source databases ==="
# Negatives come from the Metaxa2 LSU database when the download includes it;
# otherwise SILVA's LSU set is downloaded too.
if ls "$DB/LSU/blast".n* >/dev/null 2>&1; then
  bash "$ROOT/dev/sim/fetch_sources.sh" "$OUT/db" >/dev/null
  [[ -s "$OUT/metaxa2_lsu.fasta" ]] || fastacmd -d "$DB/LSU/blast" -D 1 > "$OUT/metaxa2_lsu.fasta"
  LSU="$OUT/metaxa2_lsu.fasta"; LSU_LABEL="Metaxa2 LSU"
else
  bash "$ROOT/dev/sim/fetch_sources.sh" "$OUT/db" --with-silva-lsu >/dev/null
  LSU="$OUT/db/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"; LSU_LABEL="SILVA LSU"
fi
tail -n +2 "$OUT/db/download_record.txt" | head -n 1
echo "negatives from: $LSU_LABEL"

echo "=== 2/5 sources independent of the Metaxa2 database ==="
if [[ ! -s "$OUT/sources/sources.tsv" ]]; then
  ls "$DB/SSU/blast".n* >/dev/null 2>&1 || { echo "FAIL: no BLAST database at $DB/SSU/blast" >&2; exit 1; }
  [[ -s "$OUT/metaxa2_ssu.fasta" ]] || fastacmd -d "$DB/SSU/blast" -D 1 > "$OUT/metaxa2_ssu.fasta"
  "${SIM[@]}" python3 "$ROOT/dev/sim/build_sources.py" \
    --silva-ssu "$OUT/db/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz" \
    --lsu "$LSU" --lsu-label "$LSU_LABEL" \
    --refseq "$OUT"/db/refseq/*.genomic.gbff.gz \
    --metaxa2-fasta "$OUT/metaxa2_ssu.fasta" --out "$OUT/sources" \
    --per-origin "$PER" --negatives "$NEG" --threads "$THREADS"
fi

echo "=== 3/5 simulated reads ==="
todo=()
for p in ${PROFILES//,/ }; do [[ -s "$OUT/reads/$p.info.tsv" ]] || todo+=("$p"); done
if [[ ${#todo[@]} -gt 0 ]]; then
  "${SIM[@]}" python3 "$ROOT/dev/sim/simulate_reads.py" --sources "$OUT/sources" --out "$OUT/reads" \
    --profiles "$(IFS=,; echo "${todo[*]}")" --reads-per-source "$RPS" --threads "$THREADS"
fi

: > "$OUT/SIM_REPORTS.txt"
for p in ${PROFILES//,/ }; do
  echo "=== 4/5 benchmark: $p ==="
  if [[ ! -s "$OUT/bench/$p/report.txt" ]]; then
    "$ROOT/dev/run_benchmark.sh" --db "$DB" --reads "$OUT/reads/$p.fasta" \
      --truth "$OUT/reads/$p.truth.tsv" --out "$OUT/bench/$p" --threads "$THREADS" --bin "$BIN"
  fi
  echo "=== 5/5 by novelty band: $p ==="
  python3 "$ROOT/dev/sim/score_by_band.py" --calls "$OUT/bench/$p/calls.tsv" \
    --info "$OUT/reads/$p.info.tsv" --sources "$OUT/sources/sources.tsv" | tee "$OUT/bench/$p/bands.txt"
  { echo "################ $p ################"; cat "$OUT/bench/$p/report.txt" "$OUT/bench/$p/bands.txt"; echo; } \
    >> "$OUT/SIM_REPORTS.txt"
done
cat "$OUT/bench/"*/versions.txt 2>/dev/null | head -n 7 > "$OUT/versions.txt" || true
echo; echo "Done. Combined report: $OUT/SIM_REPORTS.txt"
