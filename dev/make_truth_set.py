#!/usr/bin/env python3
"""Build a small *labelled* SSU benchmark set from Metaxa2's reference DB.

Metaxa2's SSU reference lives in a (legacy) BLAST DB under <metaxa2_db>/SSU/.
Each reference ID ends in a single-letter origin code, e.g.

    CP002869.3649570.3651094.B   ->  B = bacteria
    AY829254.1.1691.E            ->  E = eukaryota

which is the authoritative origin label (it matches the HMM codes, and is more
reliable than the taxonomy string -- a metazoan nuclear 18S is coded E, not M).

This dumps the DB with `fastacmd`, reads the origin from each ID's suffix,
reservoir-samples N per origin, and writes:

    <out>.reads.fasta   reads with clean IDs (bacteria_03, chloroplast_07, ...)
    <out>.truth.tsv     read_id <tab> true_origin

A few random "negative" sequences are appended; both tools should leave those
unclassified. Metaxa2's N (mitozoa) is folded into mitochondria, matching the
HMM assembler.

Usage:
  dev/make_truth_set.py --db <metaxa2_db> --out /tmp/bench \
      [--per-origin 10] [--negatives 10] [--seed 42] [--marker SSU]

Requires `fastacmd` (legacy BLAST, ships with the metaxa2 env) on PATH.
"""
import argparse
import glob
import os
import random
import re
import subprocess
import sys
from shutil import which

# Metaxa2 single-letter origin codes -> SSUplex origins (N folds into mito).
LETTER = {
    "A": "archaea",
    "B": "bacteria",
    "C": "chloroplast",
    "E": "eukaryota",
    "M": "mitochondria",
    "N": "mitochondria",
}
ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]


def find_blastdb(db_root, marker):
    """Find the legacy BLAST prefix under <db_root>/<marker>/ (or db_root)."""
    for d in (os.path.join(db_root, marker), db_root):
        hits = glob.glob(os.path.join(d, "*.nin"))
        if hits:
            return re.sub(r"\.nin$", "", sorted(hits)[0])
    return None


def origin_of_id(seq_id):
    """Origin from the trailing single-letter code of a reference ID."""
    last = seq_id.rsplit(".", 1)[-1].upper()
    return LETTER.get(last)


def clean_id(defline_first_token):
    """'lcl|A45315.1.1521.B' -> 'A45315.1.1521.B' (strip BLAST namespace)."""
    tok = defline_first_token.lstrip(">")
    return tok.rsplit("|", 1)[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="metaxa2_db root (contains <marker>/)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--marker", default="SSU")
    ap.add_argument("--per-origin", type=int, default=10)
    ap.add_argument("--negatives", type=int, default=10)
    ap.add_argument(
        "--fragment-length",
        type=int,
        default=0,
        help="if >0, emit a sub-window of this many bp from each reference "
        "sequence (simulates short amplicon/ONT reads) instead of the full "
        "~1.5 kb SSU. Truth labels are unchanged.",
    )
    ap.add_argument(
        "--error-rate",
        type=float,
        default=0.0,
        help="if >0, substitute each base with this probability (simulates "
        "sequencing error); applied to reference-derived reads only.",
    )
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    if not which("fastacmd"):
        sys.exit("FAIL: fastacmd not on PATH (activate the metaxa2 env).")

    prefix = find_blastdb(args.db, args.marker)
    if not prefix:
        sys.exit(f"FAIL: no BLAST DB (*.nin) under {args.db}/{args.marker} or {args.db}")
    print(f">> using BLAST DB prefix: {prefix}")

    print(">> dumping reference (fastacmd -D 1) and sampling")
    reservoir = {o: [] for o in ORIGINS}
    seen = {o: 0 for o in ORIGINS}
    n_records = 0
    unknown_letters = {}

    proc = subprocess.Popen(
        ["fastacmd", "-d", prefix, "-D", "1"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )

    cur_id = None
    cur_origin = None
    buf = []

    def flush(seq_id, origin, seqparts):
        if seq_id is None or origin is None:
            return
        seq = "".join(seqparts).replace(" ", "")
        if not seq:
            return
        seen[origin] += 1
        k = seen[origin]
        res = reservoir[origin]
        if len(res) < args.per_origin:
            res.append((seq_id, seq))
        else:
            j = rng.randint(0, k - 1)
            if j < args.per_origin:
                res[j] = (seq_id, seq)

    for line in proc.stdout:
        if line.startswith(">"):
            flush(cur_id, cur_origin, buf)
            n_records += 1
            cur_id = clean_id(line.split()[0])
            cur_origin = origin_of_id(cur_id)
            if cur_origin is None:
                suf = cur_id.rsplit(".", 1)[-1]
                unknown_letters[suf] = unknown_letters.get(suf, 0) + 1
            buf = []
        else:
            buf.append(line.strip())
    flush(cur_id, cur_origin, buf)
    err = proc.stderr.read()
    proc.wait()

    if n_records == 0:
        sys.exit(f"FAIL: fastacmd produced no records.\n{err.strip()}")

    print(f">> scanned {n_records} reference sequences")
    for o in ORIGINS:
        print(f"   {o:<13} {seen[o]} available, sampled {len(reservoir[o])}")
    if unknown_letters:
        print(f"   note: ID suffixes not mapped to an origin: {unknown_letters}")
    if all(len(v) == 0 for v in reservoir.values()):
        sys.exit("No origins recovered from ID suffixes -- paste a few deflines and we'll adapt.")

    reads_path, truth_path = args.out + ".reads.fasta", args.out + ".truth.tsv"
    counts = {o: 0 for o in ORIGINS}

    def fragment(seq):
        """Optionally take a random window of --fragment-length bp."""
        L = args.fragment_length
        if L <= 0 or len(seq) <= L:
            return seq
        start = rng.randint(0, len(seq) - L)
        return seq[start:start + L]

    def add_errors(seq):
        """Optionally substitute each base with prob --error-rate."""
        r = args.error_rate
        if r <= 0:
            return seq
        out = []
        for b in seq:
            if rng.random() < r:
                out.append(rng.choice([x for x in "ACGT" if x != b.upper()]))
            else:
                out.append(b)
        return "".join(out)

    with open(reads_path, "w") as fa, open(truth_path, "w") as tt:
        tt.write("read_id\ttrue_origin\n")
        for o in ORIGINS:
            for _, seq in reservoir[o]:
                seq = add_errors(fragment(seq))
                counts[o] += 1
                rid = f"{o}_{counts[o]:03d}"
                fa.write(f">{rid}\n")
                for i in range(0, len(seq), 80):
                    fa.write(seq[i:i + 80] + "\n")
                tt.write(f"{rid}\t{o}\n")
        for i in range(1, args.negatives + 1):
            rid = f"negative_{i:03d}"
            length = args.fragment_length if args.fragment_length > 0 else rng.randint(150, 400)
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            fa.write(f">{rid}\n{seq}\n")
            tt.write(f"{rid}\t(negative)\n")

    total = sum(counts.values()) + args.negatives
    frag = f" (fragmented to ~{args.fragment_length} bp)" if args.fragment_length > 0 else ""
    erra = f" (~{args.error_rate:.1%} substitution error)" if args.error_rate > 0 else ""
    print(f"\n>> wrote {total} reads{frag}{erra} -> {reads_path}")
    print(f">> truth labels      -> {truth_path}")


if __name__ == "__main__":
    main()
