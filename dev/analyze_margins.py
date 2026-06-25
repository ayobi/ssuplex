#!/usr/bin/env python3
"""Analyse origin-score separability from a Metaxa2 .extraction.results file.

For each read we parse Metaxa2's per-origin region tally (column 15: e.g.
``A: 10 0.0011 23.9, B: 15 2.4e-05 43.1, C: 12 1.4e-05 40.4, M: 10 0.001 45.5``)
and, grouped by the origin Metaxa2 actually called (column 3 letter code),
report the distribution of:

  * margin       = best_organelle_mean - best_prokaryote_mean
  * org_mean     = best_organelle mean per-region bit score
  * prok_mean    = best_prokaryote (A/B) mean per-region bit score
  * org_log10E   = log10 of the best organelle's E-value (more negative = stronger)

The question this answers: do reads Metaxa2 calls M/C (true organelle) have a
margin / org_mean / E-value distribution that is cleanly separated from reads it
calls B (true prokaryote, where an organelle merely "almost wins")? If the
distributions overlap, no single-threshold margin rule can separate them.

Usage:
    python3 dev/analyze_margins.py --metaxa2 bench_rice/metaxa2/run.extraction.results
"""
import argparse
import math
import statistics as st

LETTER = {"A": "archaea", "B": "bacteria", "C": "chloroplast",
          "E": "eukaryota", "M": "mito", "N": "mito"}


def parse_tally(field):
    """'A: 10 0.0011 23.9, B: 15 2.4e-05 43.1' -> {origin: (n, evalue, mean)}"""
    out = {}
    for chunk in field.split(","):
        parts = chunk.split()
        if len(parts) < 4:
            continue
        letter = parts[0].rstrip(":")
        origin = LETTER.get(letter)
        if origin is None:
            continue
        try:
            n = int(parts[1]); evalue = float(parts[2]); mean = float(parts[3])
        except ValueError:
            continue
        # fold M+N into mito, keep the stronger (higher mean)
        if origin not in out or mean > out[origin][2]:
            out[origin] = (n, evalue, mean)
    return out


def summarize(name, values):
    if not values:
        print(f"  {name:12s} (n=0)")
        return
    values = sorted(values)
    def pct(p):
        return values[min(len(values) - 1, int(p * len(values)))]
    print(f"  {name:12s} n={len(values):5d}  "
          f"min={values[0]:6.1f}  p10={pct(.10):6.1f}  median={pct(.50):6.1f}  "
          f"p90={pct(.90):6.1f}  max={values[-1]:6.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metaxa2", required=True,
                    help="Metaxa2 *.extraction.results file")
    args = ap.parse_args()

    groups = {}  # called_origin -> dict of metric -> list
    with open(args.metaxa2) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue
            called = LETTER.get(cols[2])
            if called is None:
                continue
            tally = parse_tally(cols[14])
            prok = max((tally[o][2] for o in ("archaea", "bacteria") if o in tally),
                       default=None)
            orgs = {o: tally[o] for o in ("chloroplast", "mito") if o in tally}
            if prok is None or not orgs:
                continue
            best_org_label = max(orgs, key=lambda o: orgs[o][2])
            org_n, org_e, org_mean = orgs[best_org_label]
            g = groups.setdefault(called, {"margin": [], "org_mean": [],
                                           "prok_mean": [], "org_log10E": []})
            g["margin"].append(org_mean - prok)
            g["org_mean"].append(org_mean)
            g["prok_mean"].append(prok)
            g["org_log10E"].append(math.log10(org_e) if org_e > 0 else -300.0)

    print("Distributions grouped by the origin Metaxa2 CALLED")
    print("(margin = best_organelle_mean - best_prokaryote_mean)\n")
    for called in ("bacteria", "mito", "chloroplast", "archaea", "eukaryota"):
        if called not in groups:
            continue
        print(f"=== Metaxa2 called {called} ===")
        for metric in ("margin", "org_mean", "prok_mean", "org_log10E"):
            summarize(metric, groups[called][metric])
        print()


if __name__ == "__main__":
    main()
