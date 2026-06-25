#!/usr/bin/env bash
# Assemble SSUplex-style per-origin HMM files from a Metaxa2 database.
#
# SSUplex expects one profile file per origin: bacteria.hmm, archaea.hmm,
# eukaryota.hmm, mitochondria.hmm, chloroplast.hmm. Metaxa2 ships its SSU
# profiles under <metaxa2_db>/SSU/HMMs/. This script discovers the .hmm files
# there, buckets them by origin from their path, and concatenates each bucket
# into <outdir>/<origin>.hmm.
#
# Because the exact filenames vary by Metaxa2 version, the script PRINTS what
# it found and lists anything it could not bucket — eyeball that before
# trusting the result.
#
# Usage:
#   dev/assemble_metaxa2_hmms.sh <metaxa2_db_or_HMMs_dir> <outdir> [--press]
#
# <metaxa2_db_or_HMMs_dir> may be the metaxa2_db root (the script will look for
# SSU/HMMs under it) or the HMMs directory directly.
set -euo pipefail

[[ $# -ge 2 ]] || { echo "usage: $0 <metaxa2_db_or_HMMs_dir> <outdir> [--press]"; exit 1; }
SRC="$1"; OUT="$2"; PRESS="${3:-}"

# Locate the SSU HMM directory. Prefer SSU/HMMs so that pointing at the
# metaxa2_db ROOT does not accidentally sweep in LSU/COI profiles. Search is
# shallow (maxdepth 1) and ignores macOS AppleDouble (._*) files.
HMMDIR=""
for cand in "$SRC/SSU/HMMs" "$SRC/HMMs" "$SRC"; do
  if [[ -d "$cand" ]] && [[ -n "$(find "$cand" -maxdepth 1 -name '*.hmm' ! -name '._*' -print -quit 2>/dev/null)" ]]; then
    HMMDIR="$cand"; break
  fi
done
[[ -n "$HMMDIR" ]] || { echo "FAIL: no .hmm files found under '$SRC' (tried SSU/HMMs, HMMs, then it directly)"; exit 1; }

echo ">> using HMM source: $HMMDIR"
mapfile -t ALLHMM < <(find "$HMMDIR" -maxdepth 1 -name '*.hmm' ! -name '._*' | sort)
echo ">> found ${#ALLHMM[@]} .hmm file(s):"
printf '   %s\n' "${ALLHMM[@]}"

mkdir -p "$OUT"

# origin -> regex matched against the lowercased file path.
# Metaxa2's single-letter SSU profiles: A=archaea B=bacteria C=chloroplast
# E=eukaryota M=mitochondria N=mitozoa (metazoan/animal mito SSU). N is folded
# into mitochondria for SSUplex's five-origin scheme.
declare -A PAT=(
  [bacteria]='bacteria|/b[._/]|bacterial'
  [archaea]='archaea|/a[._/]|archaeal'
  [eukaryota]='eukaryot|eukary|/e[._/]'
  [mitochondria]='mitochond|mitozoa|metazoa|/m[._/]|/n[._/]'
  [chloroplast]='chloroplast|plastid|/c[._/]'
)
ORDER=(bacteria archaea eukaryota mitochondria chloroplast)

declare -A BUCKETED
unbucketed=()
echo; echo ">> bucketing"
for f in "${ALLHMM[@]}"; do
  lc="$(echo "$f" | tr '[:upper:]' '[:lower:]')"
  hit=""
  for o in "${ORDER[@]}"; do
    if [[ "$lc" =~ ${PAT[$o]} ]]; then hit="$o"; break; fi
  done
  if [[ -n "$hit" ]]; then
    BUCKETED[$hit]="${BUCKETED[$hit]:-} $f"
  else
    unbucketed+=("$f")
  fi
done

for o in "${ORDER[@]}"; do
  files="${BUCKETED[$o]:-}"
  if [[ -z "$files" ]]; then
    echo "   $o: (no files matched — $o.hmm NOT written)"
    continue
  fi
  # shellcheck disable=SC2086
  n=$(echo $files | wc -w)
  note=""
  if [[ "$o" == "mitochondria" ]] && echo "$files" | grep -qiE '/n\.hmm'; then
    note="  (includes Metaxa2 N=mitozoa folded into mitochondria)"
  fi
  # shellcheck disable=SC2086
  cat $files > "$OUT/$o.hmm"
  echo "   $o: $n file(s) -> $OUT/$o.hmm$note"
done

if [[ ${#unbucketed[@]} -gt 0 ]]; then
  echo
  echo "!! ${#unbucketed[@]} file(s) could not be bucketed by origin — INSPECT THESE:"
  printf '   %s\n' "${unbucketed[@]}"
  echo "   (edit the PAT regexes above, or cat them into the right <origin>.hmm by hand)"
fi

if [[ "$PRESS" == "--press" ]]; then
  echo; echo ">> hmmpress-ing assembled profiles"
  for f in "$OUT"/*.hmm; do [[ -e "$f" ]] && hmmpress -f "$f" >/dev/null && echo "   pressed $f"; done
fi

echo; echo ">> done. point SSUplex at: --hmm-dir $OUT"
