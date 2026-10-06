#!/usr/bin/env bash
# Timing harness: wall time and peak memory of SSUplex and Metaxa2 on the first
# N reads of a FASTA file, for several N, with replicates. Memory is measured
# over the whole process tree (dev/peakmem.py), so helper processes such as
# nhmmer, hmmsearch and BLAST are counted.
#
# Both tools should use the same HMMER build; versions.txt records which ones
# were found. Metaxa2 runs extraction only (-x T) by default, which is the
# step SSUplex reimplements; --m2-mode full times the default pipeline too.
#
# Requires on PATH: metaxa2, nhmmer, python3 with psutil (matplotlib for the
# figure).
#
# Usage:
#   dev/run_timing.sh --reads <reads.fasta> --hmm-dir <hmms> --out <dir> \
#       [--sizes "5000 25000 100000 200000 500000 1000000"] [--reps 3] \
#       [--threads 12] [--bin path/to/ssuplex] [--m2-mode extract|full]
#
# Re-running resumes: runs already in <out>/timing.tsv are skipped.
set -euo pipefail

READS=""; HMM=""; OUT=""; SIZES="5000 25000 100000 200000 500000 1000000"
REPS=3; THREADS=12; BIN=""; M2MODE="extract"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)   READS="$2"; shift 2 ;;
    --hmm-dir) HMM="$2"; shift 2 ;;
    --out)     OUT="$2"; shift 2 ;;
    --sizes)   SIZES="$2"; shift 2 ;;
    --reps)    REPS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --bin)     BIN="$2"; shift 2 ;;
    --m2-mode) M2MODE="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$READS" && -n "$HMM" && -n "$OUT" ]] || {
  echo "usage: $0 --reads <reads.fasta> --hmm-dir <hmms> --out <dir> [--sizes ...] [--reps N] [--threads N] [--bin BIN] [--m2-mode extract|full]" >&2
  exit 1; }

[[ -s "$READS" ]] || { echo "FAIL: input FASTA not found or empty: $READS" >&2; exit 1; }
[[ -s "$HMM/bacteria.hmm" ]] || { echo "FAIL: no bacteria.hmm in $HMM" >&2; exit 1; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in metaxa2 nhmmer hmmsearch python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FAIL: '$tool' not on PATH" >&2; exit 1; }
done
python3 -c "import psutil" 2>/dev/null || { echo "FAIL: python3 cannot import psutil" >&2; exit 1; }
if [[ -z "$BIN" ]]; then
  ( cd "$ROOT" && cargo build --release --locked --quiet )
  BIN="$ROOT/target/release/ssuplex"
fi

mkdir -p "$OUT/subsets" "$OUT/runs"
{
  echo "date:      $(date)"
  echo "machine:   $(uname -srm)"
  if [[ -r /proc/cpuinfo ]]; then
    echo "cpu:       $(grep -m 1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //'), $(nproc) logical CPUs"
    echo "memory:    $(awk '/MemTotal/{printf "%.1f GB", $2/1e6}' /proc/meminfo)"
  else
    echo "cpu:       $(sysctl -n machdep.cpu.brand_string 2>/dev/null), $(sysctl -n hw.logicalcpu 2>/dev/null) logical CPUs"
    echo "memory:    $(sysctl -n hw.memsize 2>/dev/null | awk '{printf "%.1f GB", $1/1e9}')"
  fi
  echo "ssuplex:   $("$BIN" --version)"
  echo "nhmmer:    $(command -v nhmmer)  $(nhmmer -h | sed -n 2p)"
  echo "hmmsearch: $(command -v hmmsearch)  $(hmmsearch -h | sed -n 2p)"
  echo "metaxa2:   $(command -v metaxa2)  $(metaxa2 --help 2>&1 | grep -i -m 1 'version' || true)"
  echo "settings:  threads $THREADS, replicates $REPS, Metaxa2 mode $M2MODE, input $READS"
} > "$OUT/versions.txt"
cat "$OUT/versions.txt"

TSV="$OUT/timing.tsv"
[[ -s "$TSV" ]] || printf 'tool\treads\trep\tthreads\twall_s\ttree_rss_mb\tmax_rss_mb\tcpu_s\texit\n' > "$TSV"
total="$(grep -c '^>' "$READS")"

for n in $SIZES; do
  if (( n > total )); then echo "skipping $n reads: input has only $total"; continue; fi
  sub="$OUT/subsets/first_$n.fasta"
  [[ -s "$sub" ]] || awk -v n="$n" '/^>/{c++} c>n{exit} {print}' "$READS" > "$sub"
  for rep in $(seq 1 "$REPS"); do
    for tool in ssuplex metaxa2; do
      if awk -F'\t' -v t="$tool" -v n="$n" -v r="$rep" '$1==t && $2==n && $3==r {f=1} END{exit !f}' "$TSV"; then
        continue
      fi
      d="$OUT/runs/${tool}_${n}_rep${rep}"
      rm -rf "$d"; mkdir -p "$d"
      if [[ "$tool" == ssuplex ]]; then
        cmd=("$BIN" -i "$sub" -o "$d/out" --hmm-dir "$HMM" -t "$THREADS")
      else
        cmd=(metaxa2 -i "$sub" -o "$d/out" --cpu "$THREADS" -g ssu)
        [[ "$M2MODE" == extract ]] && cmd+=(-x T)
      fi
      res="$(python3 "$ROOT/dev/peakmem.py" --log "$d/log.txt" -- "${cmd[@]}")"
      printf '%s\t%s\t%s\t%s\t%s\n' "$tool" "$n" "$rep" "$THREADS" "$res" >> "$TSV"
      echo "  $tool  reads $n  rep $rep:  wall_s tree_rss_mb max_rss_mb cpu_s exit = $res"
      # Keep logs and summaries; drop large sequence and graph outputs.
      find "$d" -type f \( -name 'out*.fasta' -o -name 'out*.graph' \) -delete
    done
  done
done

python3 "$ROOT/dev/plot_timing.py" "$TSV" "$OUT/performance" \
  || echo "figure not drawn (is matplotlib installed?)"
