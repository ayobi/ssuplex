#!/usr/bin/env python3
"""Compare per-read origin calls between SSUplex and Metaxa2.

Joins the two tools' outputs by read ID and reports a confusion matrix plus
agreement statistics. Designed to be robust to Metaxa2 output-format variation:
the Metaxa2 origin label is auto-detected by scanning each row for a known
origin token rather than assuming a fixed column index.

Inputs
------
--ssuplex  SSUplex `*.extraction.tsv` (read_id, origin, score, ...)
--metaxa2     Metaxa2 `*.extraction.results` (tab-separated; id in col 1,
              origin somewhere in the row)
--reads       (optional) the input FASTA, so reads neither tool classified
              are still counted in the universe
--out         (optional) path prefix for a CSV dump of the per-read join

Interpretation note
-------------------
Metaxa2 requires >=2 conserved SSU domains to hit before it classifies a read,
which keeps its false-positive rate near zero. SSUplex v0.1 takes the single
best hit, so it is more permissive. Expect a bucket of reads SSUplex calls
to a concrete origin that Metaxa2 leaves undetected ("SSUplex-only") — those
are the calls to scrutinise, not assume correct.
"""
import argparse
import sys
from collections import Counter, defaultdict

# Canonical origins + synonyms seen across Metaxa2 versions / docs.
CANON = {
    "bacteria": "bacteria",
    "bacterial": "bacteria",
    "archaea": "archaea",
    "archaeal": "archaea",
    "eukaryota": "eukaryota",
    "eukaryote": "eukaryota",
    "eukaryotic": "eukaryota",
    "eukarya": "eukaryota",
    "nuclear": "eukaryota",
    "mitochondria": "mitochondria",
    "mitochondrial": "mitochondria",
    "mitozoa": "mitochondria",
    "metazoa": "mitochondria",
    "chloroplast": "chloroplast",
    "plastid": "chloroplast",
}
ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]

# Labels used when a tool did not assign a concrete origin.
NOT_DETECTED = "(not detected)"  # Metaxa2 produced no row for this read
UNCLASSIFIED = "unclassified"    # SSUplex's explicit unclassified bucket


def canon(token: str):
    return CANON.get(token.strip().lower())


def read_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def parse_ssuplex(path):
    """Return {read_id: origin}. Trusts the TSV's own header layout."""
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_id = header.index("read_id")
            i_origin = header.index("origin")
        except ValueError:
            sys.exit(f"unexpected SSUplex header in {path}: {header}")
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            out[fields[i_id]] = fields[i_origin].strip().lower()
    return out


# Metaxa2 .extraction.results column 3 is a single-letter origin code.
M2_LETTER = {
    "A": "archaea",
    "B": "bacteria",
    "C": "chloroplast",
    "E": "eukaryota",
    "M": "mitochondria",
    "N": "mitochondria",  # mitozoa folds into mitochondria
}


def parse_metaxa2(path):
    """Return {read_id: origin}.

    Metaxa2's `.extraction.results` is tab-separated with the origin as a
    single-letter code in column 3 (A/B/C/E/M/N). We read that positionally;
    if column 3 isn't a known letter (other Metaxa2 versions/formats), we fall
    back to scanning the row for a recognisable origin word.
    """
    out = {}
    unmatched = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            read_id = fields[0].split()[0]
            origin = None
            if len(fields) >= 3:
                origin = M2_LETTER.get(fields[2].strip().upper())
            if origin is None:  # fallback: keyword scan
                for fld in fields[1:]:
                    c = canon(fld)
                    if c:
                        origin = c
                        break
            if origin is None:
                unmatched += 1
                origin = "(unparsed)"
            out[read_id] = origin
    if unmatched:
        print(
            f"WARNING: {unmatched} Metaxa2 row(s) had no recognisable origin "
            f"(neither a column-3 letter nor a keyword) — inspect {path}.",
            file=sys.stderr,
        )
    return out


def parse_truth(path):
    """Return {read_id: true_origin}; keeps '(negative)' as-is."""
    out = {}
    with open(path) as f:
        header = f.readline()  # read_id<tab>true_origin
        _ = header
        for line in f:
            if not line.strip():
                continue
            rid, origin = line.rstrip("\n").split("\t")[:2]
            out[rid] = origin.strip().lower()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ssuplex", required=True)
    ap.add_argument("--metaxa2", required=True)
    ap.add_argument("--reads")
    ap.add_argument("--truth", help="optional truth TSV (read_id, true_origin) to score both tools against")
    ap.add_argument("--out")
    args = ap.parse_args()

    mr = parse_ssuplex(args.ssuplex)
    m2 = parse_metaxa2(args.metaxa2)

    if args.reads:
        universe = list(dict.fromkeys(read_fasta_ids(args.reads)))
    else:
        universe = list(dict.fromkeys(list(mr) + list(m2)))

    rows = []  # (read_id, m2_label, mr_label)
    for rid in universe:
        m2_label = m2.get(rid, NOT_DETECTED)
        mr_label = mr.get(rid, NOT_DETECTED)
        rows.append((rid, m2_label, mr_label))

    # Confusion matrix: rows = Metaxa2, cols = SSUplex.
    row_keys = ORIGINS + [NOT_DETECTED, "(unparsed)"]
    col_keys = ORIGINS + [UNCLASSIFIED, NOT_DETECTED]
    matrix = defaultdict(Counter)
    for _, m2l, mrl in rows:
        matrix[m2l][mrl] += 1

    def present_rows():
        return [k for k in row_keys if sum(matrix[k].values()) > 0]

    def present_cols():
        seen = set()
        for k in row_keys:
            seen.update(matrix[k].keys())
        return [c for c in col_keys if c in seen] + [
            c for c in seen if c not in col_keys
        ]

    pr, pc = present_rows(), present_cols()
    w = max([len(c) for c in pc] + [12]) + 1

    print("\nConfusion matrix  (rows = Metaxa2, cols = SSUplex)\n")
    print(" " * 15 + "".join(f"{c:>{w}}" for c in pc))
    for r in pr:
        print(f"{r:<15}" + "".join(f"{matrix[r][c]:>{w}}" for c in pc))

    # Agreement buckets.
    both_concrete = agree = disagree = 0
    mr_only = m2_only = neither = 0
    for _, m2l, mrl in rows:
        m2_concrete = m2l in ORIGINS
        mr_concrete = mrl in ORIGINS
        if m2_concrete and mr_concrete:
            both_concrete += 1
            if m2l == mrl:
                agree += 1
            else:
                disagree += 1
        elif mr_concrete and not m2_concrete:
            mr_only += 1
        elif m2_concrete and not mr_concrete:
            m2_only += 1
        else:
            neither += 1

    total = len(rows)
    pct = lambda n, d: (100.0 * n / d) if d else 0.0
    print(f"\nReads considered:            {total}")
    print(f"Both assigned a concrete origin: {both_concrete}")
    print(f"  agree:                     {agree} ({pct(agree, both_concrete):.1f}% of both-concrete)")
    print(f"  disagree:                  {disagree} ({pct(disagree, both_concrete):.1f}% of both-concrete)")
    print(f"SSUplex-only (M2 undetected): {mr_only}   <- scrutinise: possible SSUplex false positives")
    print(f"Metaxa2-only (MR unclassified):  {m2_only}   <- scrutinise: possible SSUplex misses")
    print(f"Neither classified:          {neither}")

    if args.truth:
        truth = parse_truth(args.truth)

        def score(tool_map, tool_name):
            concrete_total = concrete_correct = 0
            neg_total = neg_correct = 0
            for rid, true_o in truth.items():
                call = tool_map.get(rid, NOT_DETECTED)
                if true_o in ORIGINS:
                    concrete_total += 1
                    if call == true_o:
                        concrete_correct += 1
                elif true_o == "(negative)":
                    neg_total += 1
                    # correct = tool did NOT assign a concrete origin
                    if call not in ORIGINS:
                        neg_correct += 1
            print(f"\n{tool_name} vs truth:")
            if concrete_total:
                print(f"  origin accuracy:   {concrete_correct}/{concrete_total} "
                      f"({100.0 * concrete_correct / concrete_total:.1f}%)")
            if neg_total:
                print(f"  negatives rejected: {neg_correct}/{neg_total} "
                      f"({100.0 * neg_correct / neg_total:.1f}%)")

        print("\n" + "=" * 48)
        print("Scoring against truth labels")
        score(m2, "Metaxa2")
        score(mr, "SSUplex")

    if args.out:
        csv_path = args.out if args.out.endswith(".csv") else args.out + ".join.csv"
        with open(csv_path, "w") as f:
            f.write("read_id,metaxa2,ssuplex\n")
            for rid, m2l, mrl in rows:
                f.write(f"{rid},{m2l},{mrl}\n")
        print(f"\nwrote per-read join to {csv_path}")


if __name__ == "__main__":
    main()
