#!/usr/bin/env python3
"""Compare origin-ranking rules offline on dev/run_benchmark.sh output.

For each dataset directory, rebuilds every read's per-origin hits from
<dataset>/ssuplex/mean.regions.tsv (the exact hits SSUplex ranked) and assigns
an origin with each rule below. The mean, sum and count rules are first checked
against SSUplex's own calls, so the re-implementation is validated before the
other rules are compared.

Rules:
  sum, mean, count   as in SSUplex (mitochondria = Metaxa2's M and N sets together)
  mean_MN            mean, with M and N scored as separate origins
  normsum_MN         Metaxa2's normalised sum: summed score divided by the number
                     of profiles in the origin's set, M and N separate
  normsum            normalised sum with M and N together

Datasets with truth are scored as accuracy; datasets without truth are compared
with Metaxa2's final (BLAST-assisted) origin.

Usage: compare_rank_rules.py --hmm-src <metaxa2_db>/SSU/HMMs <dataset dir> [...]
"""
import argparse
import csv
import os
from collections import Counter, defaultdict

ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]
CODES = {"A": "archaea", "B": "bacteria", "C": "chloroplast", "E": "eukaryota", "M": "M", "N": "N"}
RULES = ["sum", "mean", "count", "mean_MN", "normsum_MN", "normsum"]


def profile_counts(hmm_src):
    counts = {}
    for code, name in CODES.items():
        with open(os.path.join(hmm_src, f"{code}.hmm")) as f:
            counts[name] = sum(1 for line in f if line.startswith("NAME"))
    counts["mitochondria"] = counts["M"] + counts["N"]
    return counts


def label(group):
    return "mitochondria" if group in ("M", "N") else group


def decide(groups, rule, nprof):
    if not groups:
        return "unclassified"
    split = rule.endswith("_MN")
    agg = defaultdict(list)
    for g, scores in groups.items():
        agg[g if split else label(g)].extend(scores)

    def key(g, s):
        n, tot = len(s), sum(s)
        if rule.startswith("mean"):
            return (tot / n, n)
        if rule.startswith("sum"):
            return (tot, n)
        if rule.startswith("count"):
            return (n, tot / n)
        return (tot / nprof[g], n)  # normsum

    best = None
    for g, s in agg.items():
        k = key(g, s)
        if best is None or k > best[0] or (k == best[0] and label(g) < label(best[1])):
            best = (k, g)
    return label(best[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm-src", required=True, help="Metaxa2 <metaxa2_db>/SSU/HMMs directory")
    ap.add_argument("datasets", nargs="+")
    a = ap.parse_args()
    nprof = profile_counts(a.hmm_src)
    print("profiles per set: " + ", ".join(f"{k} {v}" for k, v in nprof.items()))

    for ds in a.datasets:
        hits = defaultdict(lambda: defaultdict(list))
        with open(os.path.join(ds, "ssuplex", "mean.regions.tsv")) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                g = r["origin"] if r["origin"] != "mitochondria" else r["region"][:1].upper()
                hits[r["read_id"]][g].append(float(r["score"]))
        calls = list(csv.DictReader(open(os.path.join(ds, "calls.tsv")), delimiter="\t"))
        ids = [r["read_id"] for r in calls]
        truth = {r["read_id"]: r["truth"] for r in calls}
        final = {r["read_id"]: r.get("Metaxa2_final", "") for r in calls}
        hmm = {r["read_id"]: r.get("Metaxa2_HMM_stage", "") for r in calls}
        pred = {rule: {i: decide(hits.get(i, {}), rule, nprof) for i in ids} for rule in RULES}

        print("\n" + "=" * 96)
        print(f"{os.path.basename(os.path.normpath(ds))}: {len(ids)} reads")
        for rank in ("sum", "mean", "count"):
            own = {r["read_id"]: r["origin"] for r in csv.DictReader(
                open(os.path.join(ds, "ssuplex", f"{rank}.extraction.tsv")), delimiter="\t")}
            same = sum(pred[rank][i] == own.get(i, "unclassified") for i in ids)
            print(f"  check: offline {rank:<5} matches SSUplex on {same}/{len(ids)} reads")
        rows = [("Metaxa2 HMM stage", hmm), ("Metaxa2 final", final)] + [(f"SSUplex {r}", pred[r]) for r in RULES]

        if any(truth[i] for i in ids):
            ssu = [i for i in ids if truth[i] in ORIGINS]
            neg = [i for i in ids if truth[i] and truth[i] not in ORIGINS]
            classes = [c for c in ORIGINS if any(truth[i] == c for i in ssu)]
            print(f"\n  {'method':<20}{'accuracy':>15}" + "".join(f"{c[:5]:>11}" for c in classes)
                  + ("   negatives" if neg else ""))
            for name, p in rows:
                ok = sum(p[i] == truth[i] for i in ssu)
                line = f"  {name:<20}{ok:>7}/{len(ssu):<5}{100 * ok / len(ssu):5.1f}%"
                for c in classes:
                    rs = [i for i in ssu if truth[i] == c]
                    line += f"{sum(p[i] == c for i in rs):>7}/{len(rs):<3}"
                if neg:
                    line += f"   {sum(p[i] not in ORIGINS for i in neg)}/{len(neg)}"
                print(line)
        else:
            print(f"\n  {'method':<20}" + "".join(f"{c[:5]:>8}" for c in ORIGINS)
                  + "  organellar  agreement with Metaxa2 final")
            for name, p in rows:
                comp = Counter(p[i] for i in ids)
                org = comp["mitochondria"] + comp["chloroplast"]
                both = [i for i in ids if p[i] in ORIGINS and final[i] in ORIGINS]
                agree = sum(p[i] == final[i] for i in both)
                print(f"  {name:<20}" + "".join(f"{comp[c]:>8}" for c in ORIGINS)
                      + f"  {100 * org / len(ids):9.1f}%  {100 * agree / len(both):9.1f}%")


if __name__ == "__main__":
    main()
