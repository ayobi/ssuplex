#!/usr/bin/env python3
"""Score Metaxa2 and SSUplex origin calls on one dataset.

Methods compared, with identical definitions:
  Metaxa2 HMM stage        origin column of <prefix>.extraction.results
  Metaxa2 final            per-origin FASTAs from the full pipeline
                           (<prefix>.bacteria.fasta, ...; a trailing '#' marks a
                           guessed origin). Absent after an extraction-only run.
  SSUplex mean/sum/count   <ssuplex-dir>/<rank>.extraction.tsv

A read that a tool did not detect, or reported as unclassified or uncertain,
is counted as "unclassified" or "uncertain". Accuracy = correct calls / all
truth-SSU reads, so unclassified and uncertain SSU reads count as errors. Truth
negatives (e.g. LSU or non-rRNA sequences) count as correct when no concrete
origin is called.

The truth TSV has a header line and two columns, read_id and true_origin, where
true_origin is bacteria, archaea, eukaryota, mitochondria, chloroplast or
(negative).

Writes <out>/calls.tsv (one row per read: every method plus SSUplex per-origin
aggregates from <ssuplex-dir>/mean.scores.tsv) and <out>/rank_disagreements.tsv
(reads where mean, sum and count differ), and prints a report.
"""
import argparse
import glob
import os
import statistics
import sys
from collections import Counter, defaultdict

ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]
COLS = ORIGINS + ["uncertain", "unclassified"]
LETTER = {"A": "archaea", "B": "bacteria", "C": "chloroplast", "E": "eukaryota",
          "M": "mitochondria", "N": "mitochondria"}
RANKS = ["mean", "sum", "count"]


def fasta_ids(path):
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].split()[0]


def read_tsv(path):
    """Yield dict rows from a TSV with a header line."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if line.strip():
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def read_truth(path):
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            if line.strip():
                rid, origin = line.rstrip("\n").split("\t")[:2]
                out[rid] = origin.strip().lower()
    return out


def metaxa2_hmm(prefix):
    path = prefix + ".extraction.results"
    out = {}
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2].strip().upper() in LETTER:
                out[fields[0].split()[0]] = LETTER[fields[2].strip().upper()]
    return out


def metaxa2_final(prefix):
    """Return ({read_id: (origin, guessed)}, [files used]) or (None, [])."""
    out, used = {}, []
    for path in sorted(glob.glob(prefix + ".*.fasta")):
        label = os.path.basename(path)[len(os.path.basename(prefix)) + 1:-len(".fasta")].lower()
        origin = label if label in ORIGINS or label == "uncertain" else None
        if origin is None:
            continue
        used.append(os.path.basename(path))
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    out[line[1:].split()[0]] = (origin, line.rstrip().endswith("#"))
    return (out if used else None), used


def ssuplex_calls(path):
    if not os.path.exists(path):
        return None
    return {r["read_id"]: r["origin"].strip().lower() for r in read_tsv(path)}


def ssuplex_scores(path):
    out = defaultdict(dict)
    if os.path.exists(path):
        for r in read_tsv(path):
            out[r["read_id"]][r["origin"]] = (int(r["n_regions"]), float(r["mean"]), float(r["sum"]))
    return out


def ssuplex_regions(path):
    out = defaultdict(list)
    if os.path.exists(path):
        for r in read_tsv(path):
            out[r["read_id"]].append((r["origin"], r["region"], float(r["score"]),
                                      int(r["env_from"]), int(r["env_to"])))
    return out


def norm(call):
    if call in ORIGINS or call == "uncertain":
        return call
    return "unclassified"


def pct(a, b):
    return f"{100.0 * a / b:5.1f}%" if b else "   n/a"


def table(rows, cols, title, row_label="", width=13):
    print(f"\n{title}")
    print(f"  {row_label:<22}" + "".join(f"{c[:width-1]:>{width}}" for c in cols))
    for name, counts in rows:
        print(f"  {name:<22}" + "".join(f"{counts.get(c, 0):>{width}}" for c in cols))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--metaxa2-prefix", required=True)
    ap.add_argument("--ssuplex-dir", required=True)
    ap.add_argument("--truth")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    universe = list(dict.fromkeys(fasta_ids(a.reads)))
    truth = read_truth(a.truth) if a.truth else None

    methods = {}  # name -> {read_id: call}
    hmm = metaxa2_hmm(a.metaxa2_prefix)
    if hmm is not None:
        methods["Metaxa2 HMM stage"] = {r: norm(hmm.get(r)) for r in universe}
    final, used = metaxa2_final(a.metaxa2_prefix)
    guessed = {}
    if final is not None:
        methods["Metaxa2 final"] = {r: norm(final.get(r, (None,))[0]) for r in universe}
        guessed = {r: final[r][1] for r in final}
    for rank in RANKS:
        calls = ssuplex_calls(os.path.join(a.ssuplex_dir, f"{rank}.extraction.tsv"))
        if calls is not None:
            methods[f"SSUplex {rank}"] = {r: norm(calls.get(r)) for r in universe}

    print("=" * 78)
    print(f"Dataset: {a.name}   reads: {len(universe)}   truth: {'yes' if truth else 'no'}")
    print(f"Metaxa2 final-origin files: {', '.join(used) if used else 'none found (extraction-only run?)'}")
    if guessed:
        print(f"Metaxa2 final calls marked '#' (guessed origin): {sum(guessed.values())}")
    print("=" * 78)

    # ---- composition -----------------------------------------------------
    table([(m, Counter(c.values())) for m, c in methods.items()], COLS,
          "Composition (reads per call)")

    # ---- accuracy against truth -------------------------------------------
    if truth:
        print("\nAccuracy against truth (unclassified/uncertain SSU reads count as errors)")
        classes = [c for c in ORIGINS if c in truth.values()]
        negs = [r for r, t in truth.items() if t not in ORIGINS]
        for m, calls in methods.items():
            ssu = [r for r, t in truth.items() if t in ORIGINS]
            correct = sum(calls.get(r) == truth[r] for r in ssu)
            line = f"  {m:<22} accuracy {correct}/{len(ssu)} ({pct(correct, len(ssu)).strip()})"
            if negs:
                rej = sum(calls.get(r) not in ORIGINS for r in negs)
                line += f"   negatives rejected {rej}/{len(negs)}"
            print(line)
            recalls = []
            for c in classes:
                rs = [r for r in ssu if truth[r] == c]
                recalls.append(f"{c[:4]} {sum(calls.get(r) == c for r in rs)}/{len(rs)}")
            print("      per-class recall: " + "   ".join(recalls))
        for m, calls in methods.items():
            rows = []
            for c in classes + (["(negative)"] if negs else []):
                rs = [r for r, t in truth.items() if (t == c if c in ORIGINS else t not in ORIGINS)]
                rows.append((f"truth {c}", Counter(calls.get(r) for r in rs)))
            table(rows, COLS, f"Confusion matrix: {m} (rows = truth, cols = call)")

    # ---- agreement with Metaxa2 -------------------------------------------
    for ref in ("Metaxa2 HMM stage", "Metaxa2 final"):
        if ref not in methods:
            continue
        print(f"\nAgreement with {ref} (among reads both assign a concrete origin)")
        for rank in RANKS:
            m = f"SSUplex {rank}"
            if m not in methods:
                continue
            both = [r for r in universe if methods[ref][r] in ORIGINS and methods[m][r] in ORIGINS]
            agree = sum(methods[ref][r] == methods[m][r] for r in both)
            print(f"  {m:<15} {agree}/{len(both)} ({pct(agree, len(both)).strip()})")
    for ref, m in (("Metaxa2 final", "SSUplex mean"), ("Metaxa2 final", "SSUplex sum")):
        if ref in methods and m in methods:
            rows = [(f"{ref.split()[1]} {c}", Counter(methods[m][r] for r in universe if methods[ref][r] == c))
                    for c in COLS]
            table(rows, COLS, f"Cross-tab: rows = {ref}, cols = {m}")

    # ---- ranking-statistic disagreements (Reviewer 2, a2) -------------------
    scores = ssuplex_scores(os.path.join(a.ssuplex_dir, "mean.scores.tsv"))
    regions = ssuplex_regions(os.path.join(a.ssuplex_dir, "mean.regions.tsv"))
    ranks_present = [r for r in RANKS if f"SSUplex {r}" in methods]
    dis = []
    if len(ranks_present) == 3:
        dis = [r for r in universe if len({methods[f"SSUplex {k}"][r] for k in RANKS}) > 1]
        print(f"\nReads where mean/sum/count disagree: {len(dis)} of {len(universe)}")
        pat = Counter(tuple(methods[f"SSUplex {k}"][r] for k in RANKS) for r in dis)
        for (mn, sm, ct), n in pat.most_common(10):
            line = f"  mean={mn:<13} sum={sm:<13} count={ct:<13} {n:>6}"
            if truth:
                tc = Counter(truth.get(r, "?") for r in dis
                             if (methods["SSUplex mean"][r], methods["SSUplex sum"][r],
                                 methods["SSUplex count"][r]) == (mn, sm, ct))
                line += "   truth: " + ", ".join(f"{k} {v}" for k, v in tc.most_common(3))
            print(line)
        if scores and dis:
            print("  Median per-origin aggregates on these reads (n regions / mean / sum):")
            for o in ORIGINS:
                vals = [scores[r][o] for r in dis if o in scores.get(r, {})]
                if vals:
                    print(f"    {o:<13} reads {len(vals):>6}   n {statistics.median(v[0] for v in vals):>5.1f}"
                          f"   mean {statistics.median(v[1] for v in vals):>7.1f}"
                          f"   sum {statistics.median(v[2] for v in vals):>8.1f}")

    # ---- mitochondria profile fold check (M + N concatenated) ---------------
    if regions:
        mito_names = Counter(reg[1][:1] for hits in regions.values() for reg in hits if reg[0] == "mitochondria")
        max_n = {o: max((v[o][0] for v in scores.values() if o in v), default=0) for o in ORIGINS}
        both_mn = overlap = 0
        for hits in regions.values():
            m = [h for h in hits if h[0] == "mitochondria" and h[1][:1] == "M"]
            n = [h for h in hits if h[0] == "mitochondria" and h[1][:1] == "N"]
            if m and n:
                both_mn += 1
                if any(min(x[4], y[4]) - max(x[3], y[3]) > 0.5 * min(x[4] - x[3], y[4] - y[3])
                       for x in m for y in n):
                    overlap += 1
        print("\nMitochondrial profile fold check (mitochondria.hmm = M.hmm + N.hmm)")
        print(f"  first letter of mitochondrial region names: {dict(mito_names)}")
        print(f"  reads with both M* and N* mitochondrial hits: {both_mn}; "
              f"with an M and N hit overlapping the same stretch: {overlap}")
        print("  max regions matched per origin: " + ", ".join(f"{o} {max_n[o]}" for o in ORIGINS))

    # ---- per-read tables ----------------------------------------------------
    os.makedirs(a.out, exist_ok=True)
    hdr = ["read_id", "truth"] + [m.replace(" ", "_") for m in methods] + ["metaxa2_final_guessed"]
    hdr += [f"{o}_{k}" for o in ORIGINS for k in ("n", "mean", "sum")]
    for fname, ids in (("calls.tsv", universe), ("rank_disagreements.tsv", dis)):
        with open(os.path.join(a.out, fname), "w") as f:
            f.write("\t".join(hdr) + "\n")
            for r in ids:
                row = [r, truth.get(r, "") if truth else ""] + [methods[m][r] for m in methods]
                row.append("yes" if guessed.get(r) else "")
                for o in ORIGINS:
                    v = scores.get(r, {}).get(o)
                    row += [str(v[0]), f"{v[1]:.2f}", f"{v[2]:.2f}"] if v else ["", "", ""]
                f.write("\t".join(row) + "\n")
    print(f"\nWrote {os.path.join(a.out, 'calls.tsv')} and rank_disagreements.tsv ({len(dis)} reads)")


if __name__ == "__main__":
    sys.exit(main())
