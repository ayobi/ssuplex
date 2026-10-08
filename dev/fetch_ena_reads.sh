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
echo "(progress is printed every 100,000 reads; large requests take a while)" >&2
# head stops the download on purpose once N reads are in, so curl's write
# error at that point is expected and silenced; the read count is checked below.
set +o pipefail
curl -fs "https://$URL" | gunzip -c 2>/dev/null | head -n $((N * 4)) \
  | awk 'NR%4==1{split($0,a," "); print ">" substr(a[1],2); if ((NR+3)/4 % 100000 == 0) printf "  %d reads\n", (NR+3)/4 > "/dev/stderr"}
         NR%4==2{print}' > "$OUT.tmp"
set -o pipefail
got="$(grep -c '^>' "$OUT.tmp" || true)"
[ "$got" -gt 0 ] || { rm -f "$OUT.tmp"; echo "FAIL: no reads downloaded for $ACC" >&2; exit 1; }
mv "$OUT.tmp" "$OUT"
if [ "$got" -lt "$N" ]; then
  echo "wrote $OUT ($got reads; the run has fewer than $N, or the download stopped early)" >&2
else
  echo "wrote $OUT ($got reads)" >&2
fi
