#!/usr/bin/env bash
# Download the source databases for the simulated benchmark into <out_dir>:
# SILVA 138.2 SSU Ref NR99 and the NCBI RefSeq mitochondrion and plastid genome
# releases (GenBank flat files), plus SILVA 138.2 LSU Ref NR99 with --with-silva-lsu
# (only needed when no other LSU source is available for negatives).
# Re-running skips completed downloads.
#
# Usage: dev/sim/fetch_sources.sh <out_dir> [--with-silva-lsu]
set -euo pipefail
OUT="${1:?usage: $0 <out_dir> [--with-silva-lsu]}"
WITH_LSU="${2:-}"
mkdir -p "$OUT/refseq"

get() {  # get <url> <dest>
  [ -s "$2" ] && return 0
  echo "downloading $1" >&2
  curl -fL --retry 3 -o "$2.part" "$1"
  mv "$2.part" "$2"
}

SILVA=https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports
get "$SILVA/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz" "$OUT/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
if [ "$WITH_LSU" = "--with-silva-lsu" ]; then
  get "$SILVA/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz" "$OUT/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
fi

REFSEQ=https://ftp.ncbi.nlm.nih.gov/refseq/release
get "$REFSEQ/RELEASE_NUMBER" "$OUT/refseq/RELEASE_NUMBER"
for kind in mitochondrion plastid; do
  files="$(curl -fsSL "$REFSEQ/$kind/" | grep -o "$kind\.[0-9]*\.genomic\.gbff\.gz" | sort -u)"
  [ -n "$files" ] || { echo "FAIL: no $kind files listed at $REFSEQ/$kind/" >&2; exit 1; }
  for f in $files; do
    get "$REFSEQ/$kind/$f" "$OUT/refseq/$f"
  done
done

{
  echo "downloaded: $(date)"
  echo "RefSeq release: $(cat "$OUT/refseq/RELEASE_NUMBER")"
  (cd "$OUT" && ls -l ./*.gz refseq/*.gz)
} > "$OUT/download_record.txt"
cat "$OUT/download_record.txt"
