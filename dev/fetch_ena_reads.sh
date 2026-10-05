#!/usr/bin/env bash
# Write the first N reads of an ENA/SRA run as FASTA (IDs truncated at the
# first space). The FASTQ location is looked up through the ENA portal API, so
# no SRA toolkit is needed; only the first N reads are downloaded.
#
# Usage: dev/fetch_ena_reads.sh <run accession> <N> <out.fasta>
#   e.g. dev/fetch_ena_reads.sh SRR10391201 200000 zymo_200k.fasta
set -euo pipefail

[ "$#" -eq 3 ] || { echo "usage: $0 <run accession> <N> <out.fasta>" >&2; exit 1; }
ACC="$1"; N="$2"; OUT="$3"

URL="$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp&format=tsv" \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="fastq_ftp") c=i; next} NR==2 && c{split($c,a,";"); print a[1]}')"
[ -n "$URL" ] || { echo "FAIL: no FASTQ found for $ACC in ENA" >&2; exit 1; }

echo "streaming the first $N reads of https://$URL" >&2
set +o pipefail   # head closes the stream early on purpose
curl -fsS "https://$URL" | gunzip -c 2>/dev/null | head -n $((N * 4)) \
  | awk 'NR%4==1{split($0,a," "); print ">" substr(a[1],2)} NR%4==2{print}' > "$OUT.tmp"
set -o pipefail
mv "$OUT.tmp" "$OUT"
echo "wrote $OUT ($(grep -c '^>' "$OUT") reads)" >&2
