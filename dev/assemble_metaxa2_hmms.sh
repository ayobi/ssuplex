#!/usr/bin/env bash
# Assemble SSUplex's five per-origin profile files from a Metaxa2 database.
#
# Usage: dev/assemble_metaxa2_hmms.sh <metaxa2_db> <out_dir>
#
# Metaxa2 stores its SSU profiles as one file per origin code in
# <metaxa2_db>/SSU/HMMs/: A archaea, B bacteria, C chloroplast, E eukaryota,
# M mitochondria and N metazoan mitochondria. SSUplex scores M and N together
# as a single mitochondrial origin. Works with the bash that ships with macOS.
set -euo pipefail

[ "$#" -eq 2 ] || { echo "usage: $0 <metaxa2_db> <out_dir>" >&2; exit 1; }
SRC="$1/SSU/HMMs"
OUT="$2"
for code in A B C E M N; do
  [ -s "$SRC/$code.hmm" ] || { echo "FAIL: $SRC/$code.hmm not found (is $1 a Metaxa2 database directory?)" >&2; exit 1; }
done

mkdir -p "$OUT"
cat "$SRC/B.hmm"              > "$OUT/bacteria.hmm"
cat "$SRC/A.hmm"              > "$OUT/archaea.hmm"
cat "$SRC/E.hmm"              > "$OUT/eukaryota.hmm"
cat "$SRC/M.hmm" "$SRC/N.hmm" > "$OUT/mitochondria.hmm"
cat "$SRC/C.hmm"              > "$OUT/chloroplast.hmm"

echo "Assembled profiles in $OUT:"
for o in bacteria archaea eukaryota mitochondria chloroplast; do
  printf '  %-13s %3d profiles\n' "$o" "$(grep -c '^NAME' "$OUT/$o.hmm")"
done
