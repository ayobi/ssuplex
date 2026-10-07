#!/usr/bin/env python3
"""Break down accuracy on simulated reads by novelty band and source subgroup.

Joins the per-read calls from dev/run_benchmark.sh (calls.tsv) with the simulation
record (<profile>.info.tsv) and the sources table (sources.tsv), then reports for
every method:
  - recall per origin within each novelty band (identity of the source to its
    closest Metaxa2 database sequence: 99% or more, 97 to 99%, below 97%)
  - recall for metazoan and non-metazoan mitochondria
  - negatives rejected
  - expected results for two sample compositions, computed from the per-origin
    call rates (each read is classified independently, so per-origin rates do not
    depend on the mix): overall accuracy, and the organellar share the method would
    report against the true share

Usage: score_by_band.py --calls <benchmark>/calls.tsv --info <profile>.info.tsv \
                        --sources <sources>/sources.tsv
"""
import argparse
import csv
from collections import Counter, defaultdict

ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]
BANDS = ["99-100", "97-99", "<97"]
MIXES = {
    "bacteria-dominated": {"bacteria": 0.90, "archaea": 0.02, "eukaryota": 0.02,
                           "mitochondria": 0.03, "chloroplast": 0.03},
    "plant-root-like": {"bacteria": 0.55, "archaea": 0.01, "eukaryota": 0.02,
                        "mitochondria": 0.07, "chloroplast": 0.35},
}
METHODS = ["Metaxa2_HMM_stage", "Metaxa2_final", "SSUplex_sum", "SSUplex_mean", "SSUplex_count"]


def table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pct(a, b):
    return f"{100 * a / b:5.1f}" if b else "    -"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--info", required=True)
    ap.add_argument("--sources", required=True)
    a = ap.parse_args()

    src = {r["source_id"]: r for r in table(a.sources)}
    read_src = {r["read_id"]: r["source_id"] for r in table(a.info)}
    calls = table(a.calls)
    methods = [m for m in METHODS if m in calls[0]]

    reads = []
    for r in calls:
        s = src.get(read_src.get(r["read_id"], ""), None)
        reads.append({"truth": r["truth"], "band": s["band"] if s else "",
                      "subgroup": s["subgroup"] if s else "", **{m: r[m] for m in methods}})
    ssu = [r for r in reads if r["truth"] in ORIGINS]
    neg = [r for r in reads if r["truth"] not in ORIGINS]
    short = {m: m.replace("Metaxa2_HMM_stage", "M2 HMM").replace("Metaxa2_final", "M2 final")
             .replace("SSUplex_", "SSUplex ") for m in methods}

    print(f"reads: {len(ssu)} SSU, {len(neg)} negatives")
    head = f"  {'origin':<13}{'band':<8}{'reads':>7}" + "".join(f"{short[m]:>14}" for m in methods)
    print("\nRecall (%) by origin and novelty band")
    print(head)
    for o in ORIGINS:
        for b in BANDS + ["all"]:
            rs = [r for r in ssu if r["truth"] == o and (b == "all" or r["band"] == b)]
            if not rs:
                continue
            print(f"  {o:<13}{b:<8}{len(rs):>7}" + "".join(
                f"{pct(sum(r[m] == o for r in rs), len(rs)):>14}" for m in methods))

    print("\nMitochondria by subgroup, recall (%)")
    for g in ("metazoan", "non-metazoan"):
        rs = [r for r in ssu if r["truth"] == "mitochondria" and r["subgroup"] == g]
        if rs:
            print(f"  {g:<21}{len(rs):>7}" + "".join(
                f"{pct(sum(r[m] == 'mitochondria' for r in rs), len(rs)):>14}" for m in methods))

    if neg:
        print("\nNegatives rejected (%)")
        print(f"  {'':<21}{len(neg):>7}" + "".join(
            f"{pct(sum(r[m] not in ORIGINS for r in neg), len(neg)):>14}" for m in methods))

    # per-origin call rates
    rate = {m: {} for m in methods}
    for m in methods:
        for o in ORIGINS:
            rs = [r for r in ssu if r["truth"] == o]
            n = len(rs)
            rate[m][o] = {c: v / n for c, v in Counter(r[m] for r in rs).items()} if n else {}
    print("\nExpected results for two sample compositions")
    for mix, w in MIXES.items():
        true_org = w["mitochondria"] + w["chloroplast"]
        print(f"  {mix} (true organellar share {100 * true_org:.0f}%)")
        for m in methods:
            if any(not rate[m][o] for o in ORIGINS if w[o] > 0):
                continue
            acc = sum(w[o] * rate[m][o].get(o, 0.0) for o in ORIGINS)
            org = sum(w[o] * (rate[m][o].get("mitochondria", 0.0) + rate[m][o].get("chloroplast", 0.0))
                      for o in ORIGINS)
            print(f"    {short[m]:<15} accuracy {100 * acc:5.1f}%   reported organellar share {100 * org:5.1f}%")


if __name__ == "__main__":
    main()
