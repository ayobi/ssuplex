#!/usr/bin/env python3
"""Margin separability on SSUplex's OWN per-origin scores.

Companion to analyze_margins.py (which used Metaxa2's tallies). This reads the
`--debug-scores` table SSUplex writes ({prefix}.scores.tsv) and groups reads
by what Metaxa2 called (from concordance.join.csv), so we can see whether
SSUplex's own mean per-region scores separate true-organelle from
true-prokaryote reads at margin = 0, the way Metaxa2's tallies did.

If they DO: the fix lives in ranking (rank by mean ~= margin-at-zero), and the
Zymo failure was about which regions got counted, not the metric.
If they DON'T: SSUplex's per-region scoring on noisy reads is itself the
problem (region inclusion / E-value gating), which is the real lever.

Usage:
    python3 dev/analyze_margins_mr.py \
        --scores /tmp/rice_dbg/sample.scores.tsv \
        --join   bench_rice/concordance.join.csv
"""
import argparse

PROK = {"bacteria", "archaea"}
ORG = {"chloroplast", "mito", "mitochondria"}


def load_scores(path):
    """read_id -> {origin: mean}"""
    reads = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= idx["mean"]:
                continue
            reads.setdefault(c[idx["read_id"]], {})[c[idx["origin"]]] = float(c[idx["mean"]])
    return reads


def load_join(path):
    """read_id -> metaxa2 call. Accepts csv or tsv with a header."""
    out = {}
    with open(path) as f:
        first = f.readline()
        sep = "," if first.count(",") >= first.count("\t") else "\t"
        header = first.rstrip("\n").split(sep)
        # expect columns read_id, metaxa2, ssuplex (any order, by name)
        try:
            ri, mi = header.index("read_id"), header.index("metaxa2")
        except ValueError:
            ri, mi = 0, 1
        for line in f:
            c = line.rstrip("\n").split(sep)
            if len(c) > max(ri, mi):
                out[c[ri]] = c[mi]
    return out


def summarize(name, values):
    if not values:
        print(f"    {name:10s} (n=0)")
        return
    values = sorted(values)
    def pct(p):
        return values[min(len(values) - 1, int(p * len(values)))]
    neg = sum(1 for v in values if v < 0)
    print(f"    {name:10s} n={len(values):5d}  min={values[0]:6.1f}  "
          f"p10={pct(.10):6.1f}  median={pct(.50):6.1f}  p90={pct(.90):6.1f}  "
          f"max={values[-1]:6.1f}   ({100*neg/len(values):.0f}% negative)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="SSUplex {prefix}.scores.tsv")
    ap.add_argument("--join", required=True, help="concordance.join.csv (for Metaxa2 grouping)")
    args = ap.parse_args()

    scores = load_scores(args.scores)
    calls = load_join(args.join)

    groups = {}
    for rid, origins in scores.items():
        called = calls.get(rid)
        if called is None:
            continue
        prok = max((m for o, m in origins.items() if o in PROK), default=None)
        org = max((m for o, m in origins.items() if o in ORG), default=None)
        if prok is None or org is None:
            continue
        groups.setdefault(called, []).append(org - prok)

    print("SSUplex's own margin (best_organelle_mean - best_prokaryote_mean)")
    print("grouped by what Metaxa2 called:\n")
    for called in ("bacteria", "mitochondria", "mito", "chloroplast", "archaea", "eukaryota"):
        if called in groups:
            print(f"  === Metaxa2 called {called} ===")
            summarize("margin", groups[called])
            print()


if __name__ == "__main__":
    main()
