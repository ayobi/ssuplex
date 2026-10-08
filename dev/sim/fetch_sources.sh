#!/usr/bin/env bash
# Download the source databases for the simulated benchmark into <out_dir>:
# SILVA 138.2 SSU Ref NR99 always; the NCBI RefSeq mitochondrion and plastid genome
# releases with --with-refseq; SILVA 138.2 LSU Ref NR99 with --with-silva-lsu.
# Re-running skips completed downloads.
#
# Usage: dev/sim/fetch_sources.sh <out_dir> [--with-refseq] [--with-silva-lsu]
set -euo pipefail
OUT="${1:?usage: $0 <out_dir> [--with-refseq] [--with-silva-lsu]}"
shift
WITH_LSU=""; WITH_REFSEQ=""
for flag in "$@"; do
  case "$flag" in
    --with-silva-lsu) WITH_LSU=1 ;;
    --with-refseq)    WITH_REFSEQ=1 ;;
    *) echo "unknown option: $flag" >&2; exit 1 ;;
  esac
done
mkdir -p "$OUT/refseq"

get() {  # get <url> <dest>
  [ -s "$2" ] && return 0
  echo "downloading $1" >&2
  curl -fL --retry 3 -o "$2.part" "$1"
  mv "$2.part" "$2"
}

SILVA=https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports
get "$SILVA/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz" "$OUT/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
if [ -n "$WITH_LSU" ]; then
  get "$SILVA/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz" "$OUT/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
fi

REFSEQ=https://ftp.ncbi.nlm.nih.gov/refseq/release
if [ -n "$WITH_REFSEQ" ]; then get "$REFSEQ/RELEASE_NUMBER" "$OUT/refseq/RELEASE_NUMBER"; fi
for kind in $([ -n "$WITH_REFSEQ" ] && echo mitochondrion plastid); do
  files="$(curl -fsSL "$REFSEQ/$kind/" | grep -o "$kind\.[0-9]*\.genomic\.gbff\.gz" | sort -u)"
  [ -n "$files" ] || { echo "FAIL: no $kind files listed at $REFSEQ/$kind/" >&2; exit 1; }
  for f in $files; do
    get "$REFSEQ/$kind/$f" "$OUT/refseq/$f"
  done
done

{
  echo "downloaded: $(date)"
  if [ -n "$WITH_REFSEQ" ]; then echo "RefSeq release: $(cat "$OUT/refseq/RELEASE_NUMBER")"; fi
  (cd "$OUT" && ls -l ./*.gz refseq/*.gz 2>/dev/null) || true
} > "$OUT/download_record.txt"
cat "$OUT/download_record.txt"
