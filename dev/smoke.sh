#!/usr/bin/env bash
# End-to-end smoke test for SSUplex.
#
# Generates synthetic HMM profiles + reads, runs the binary, and asserts the
# pipeline routes reads correctly AND accounts for every input read (the
# no-hit reads must land in unclassified.fasta and the counts must add up).
#
# Requires: cargo, python3, and HMMER (hmmbuild + hmmsearch) on PATH.
#
# Usage:
#   dev/smoke.sh                 # builds release binary, then tests
#   dev/smoke.sh path/to/binary  # tests an already-built binary
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# --- locate binary ---------------------------------------------------------
if [[ $# -ge 1 ]]; then
  BIN="$1"
else
  echo ">> building release binary"
  cargo build --release
  BIN="target/release/ssuplex"
fi
[[ -x "$BIN" ]] || { echo "FAIL: binary not found/executable: $BIN"; exit 1; }

for tool in hmmbuild hmmsearch python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FAIL: '$tool' not on PATH"; exit 1; }
done

# --- generate fixtures -----------------------------------------------------
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
python3 dev/gen_fixtures.py "$WORK"

mkdir -p "$WORK/hmms" "$WORK/out"
hmmbuild --dna "$WORK/hmms/bacteria.hmm" "$WORK/bact.afa" >/dev/null
hmmbuild --dna "$WORK/hmms/archaea.hmm"  "$WORK/arch.afa" >/dev/null

# --- run -------------------------------------------------------------------
PREFIX="$WORK/out/sample01"
echo ">> running ssuplex"
"$BIN" -i "$WORK/reads.fasta" -o "$PREFIX" --hmm-dir "$WORK/hmms" -t 4 -E 1e-3

# --- assertions ------------------------------------------------------------
fails=0
check() { # check <description> <test-command...>
  local desc="$1"; shift
  if "$@"; then
    echo "  PASS: $desc"
  else
    echo "  FAIL: $desc"
    fails=$((fails + 1))
  fi
}

nrec() { grep -c '^>' "$1" 2>/dev/null || echo 0; }

echo ">> checking outputs"
check "bacteria.fasta has 2 records"       test "$(nrec "$PREFIX.bacteria.fasta")" -eq 2
check "archaea.fasta has 2 records"        test "$(nrec "$PREFIX.archaea.fasta")"  -eq 2
check "unclassified.fasta has 2 records"   test "$(nrec "$PREFIX.unclassified.fasta")" -eq 2
check "read_junk_1 captured (not dropped)" grep -q '>read_junk_1'  "$PREFIX.unclassified.fasta"
check "read_junk_2 captured (not dropped)" grep -q '>read_junk_2'  "$PREFIX.unclassified.fasta"
check "summary: classified = 4"            grep -Eq '^classified reads: +4$'   "$PREFIX.summary.txt"
check "summary: unclassified = 2"          grep -Eq '^unclassified reads: +2$' "$PREFIX.summary.txt"
check "TSV records junk reads"             grep -q 'read_junk_1.*unclassified' "$PREFIX.extraction.tsv"
# Envelope coordinates are predictable from the known flank lengths.
check "read_bact_1 clipped to coords 61-150" grep -q 'read_bact_1.*coords=61-150' "$PREFIX.bacteria.fasta"
check "read_arch_1 clipped to coords 101-190" grep -q 'read_arch_1.*coords=101-190' "$PREFIX.archaea.fasta"

echo ">> checking overwrite protection"
if "$BIN" -i "$WORK/reads.fasta" -o "$PREFIX" --hmm-dir "$WORK/hmms" -t 4 -E 1e-3 >/dev/null 2>&1; then
  echo "  FAIL: second run without --force overwrote existing outputs"; fails=$((fails + 1))
else
  echo "  PASS: second run without --force refused to overwrite"
fi
echo ">stale" > "$PREFIX.chloroplast.fasta"   # file a previous run might have left
"$BIN" -i "$WORK/reads.fasta" -o "$PREFIX" --hmm-dir "$WORK/hmms" -t 4 -E 1e-3 --force >/dev/null 2>&1
check "--force run succeeded"                test "$(nrec "$PREFIX.bacteria.fasta")" -eq 2
check "--force removed stale per-origin file" test ! -e "$PREFIX.chloroplast.fasta"

echo
if [[ "$fails" -eq 0 ]]; then
  echo "ALL CHECKS PASSED"
else
  echo "$fails CHECK(S) FAILED"
  exit 1
fi
