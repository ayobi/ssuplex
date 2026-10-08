#!/usr/bin/env python3
"""Build the manuscript's supplementary tables as LaTeX from benchmark outputs.

Reads the per-read tables (calls.tsv) written by dev/run_benchmark.sh, either as
run by dev/run_real_data.sh (--accuracy) or by dev/sim/run_simulation.sh (--sim),
and writes:

  tableS1.tex  confusion matrices for every labelled dataset and method, and the
               rice-root comparison with Metaxa2's BLAST-assisted calls
  tableS2.tex  reads on which mean, sum and count disagree: each pattern, its
               true origin (or Metaxa2's final call for rice), and the median
               region count, mean and summed bit score of the origins chosen by
               the mean and by the sum
  tableS3.tex  simulated benchmark by read accuracy, origin and novelty band
               (only with --sim)

Usage: supplementary_tables.py --accuracy <run_real_data.sh out dir>
                               [--sim <run_simulation.sh out dir>] --out <dir>
"""
import argparse
import csv
import os
import statistics
from collections import Counter

ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]
CALLS = ORIGINS + ["uncertain", "unclassified"]
METHODS = [("Metaxa2_HMM_stage", "Metaxa2 HMM"), ("Metaxa2_final", "Metaxa2 BLAST"),
           ("SSUplex_sum", "SSUplex sum"), ("SSUplex_mean", "SSUplex mean"),
           ("SSUplex_count", "SSUplex count")]
SHORT = {"bacteria": "Bact.", "archaea": "Arch.", "eukaryota": "Euk.", "mitochondria": "Mito.",
         "chloroplast": "Chlo.", "uncertain": "Unc.", "unclassified": "Uncl.", "(negative)": "Negative"}
LABELLED = [("testfasta", "Metaxa2 test file"), ("ref", "Metaxa2 reference set"),
            ("zymo", "ZymoBIOMICS Nanopore, first 5{,}000 reads")]


def load(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    return f"{x:,}".replace(",", "{,}")


def confusion(rows, label):
    present = [c for c in ORIGINS if any(r["truth"] == c for r in rows)]
    if any(r["truth"] and r["truth"] not in ORIGINS for r in rows):
        present.append("(negative)")
    out = [r"\begin{table}[htbp]", r"\centering", r"\small",
           rf"\caption*{{\textbf{{{label}.}} Rows: true origin; columns: call. Unc., uncertain "
           r"(Metaxa2 only); Uncl., undetected or unclassified.}",
           r"\begin{tabular}{llrrrrrrr}", r"\toprule",
           "Method & True origin & " + " & ".join(SHORT[c] for c in CALLS) + r" \\", r"\midrule"]
    for i, (col, name) in enumerate(METHODS):
        for j, t in enumerate(present):
            rs = [r for r in rows if (r["truth"] == t if t in ORIGINS else r["truth"] not in ORIGINS)]
            cnt = Counter(r[col] for r in rs)
            cells = [num(cnt.get(c, 0)) if cnt.get(c, 0) else "" for c in CALLS]
            out.append(f"{name if j == 0 else ''} & {SHORT.get(t, t)} & " + " & ".join(cells) + r" \\")
        if i < len(METHODS) - 1:
            out.append(r"\midrule")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(out)


def rice_table(rows):
    out = [r"\begin{table}[htbp]", r"\centering", r"\small",
           r"\caption*{\textbf{Rice-root Nanopore 16S sample (SRR25243163), first 5{,}000 reads.} "
           r"No per-read truth: reads assigned to each origin, the organellar share, and agreement "
           r"with Metaxa2's BLAST-assisted calls among reads both methods assign to an SSU origin.}",
           r"\begin{tabular}{lrrrrrrr}", r"\toprule",
           r"Method & Bact. & Arch. & Euk. & Mito. & Chlo. & Organellar (\%) & Agreement (\%) \\",
           r"\midrule"]
    for col, name in METHODS:
        cnt = Counter(r[col] for r in rows)
        org = 100 * (cnt["mitochondria"] + cnt["chloroplast"]) / len(rows)
        both = [r for r in rows if r[col] in ORIGINS and r["Metaxa2_final"] in ORIGINS]
        agree = "" if col == "Metaxa2_final" else f"{100 * sum(r[col] == r['Metaxa2_final'] for r in both) / len(both):.1f}"
        out.append(f"{name} & " + " & ".join(num(cnt[c]) for c in ORIGINS) + f" & {org:.1f} & {agree}" + r" \\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(out)


def median(vals):
    return statistics.median(vals) if vals else float("nan")


def agg(r, origin):
    n = r.get(f"{origin}_n", "")
    return (int(n), float(r[f"{origin}_mean"]), float(r[f"{origin}_sum"])) if n else None


def disagreement_table(name, rows, ref_col, min_n=2):
    dis = [r for r in rows if len({r["SSUplex_mean"], r["SSUplex_sum"], r["SSUplex_count"]}) > 1]
    pats = Counter((r["SSUplex_mean"], r["SSUplex_sum"], r["SSUplex_count"]) for r in dis)
    out = [rf"\multicolumn{{11}}{{l}}{{\textbf{{{name}}}: {num(len(dis))} of {num(len(rows))} reads}} \\"]
    other = 0
    for (mn, sm, ct), n in pats.most_common():
        if n < min_n:
            other += n
            continue
        rs = [r for r in dis if (r["SSUplex_mean"], r["SSUplex_sum"], r["SSUplex_count"]) == (mn, sm, ct)]
        ref = Counter(r[ref_col] for r in rs).most_common(2)
        ref_txt = ", ".join(f"{SHORT.get(k, k)} {num(v)}" for k, v in ref)
        cells = []
        for chosen in (mn, sm):
            a = [agg(r, chosen) for r in rs if chosen in ORIGINS]
            a = [x for x in a if x]
            cells += [f"{median([x[0] for x in a]):.0f}", f"{median([x[1] for x in a]):.1f}",
                      f"{median([x[2] for x in a]):.0f}"] if a else ["", "", ""]
        out.append(f"{SHORT[mn]} & {SHORT[sm]} & {SHORT[ct]} & {num(n)} & {ref_txt} & " + " & ".join(cells) + r" \\")
    if other:
        out.append(rf"\multicolumn{{3}}{{l}}{{other patterns}} & {num(other)} & & & & & & & \\")
    return out


MIXES = [("bacteria-dominated", {"bacteria": 0.90, "archaea": 0.02, "eukaryota": 0.02,
                                 "mitochondria": 0.03, "chloroplast": 0.03}),
         ("plant-root-like", {"bacteria": 0.55, "archaea": 0.01, "eukaryota": 0.02,
                              "mitochondria": 0.07, "chloroplast": 0.35})]
PROFILES = [("ont_r9", "92\\%"), ("ont_r10", "98.5\\%"), ("hifi_like", "99.9\\%")]
BANDS = [("99-100", "99--100"), ("97-99", "97--99"), ("<97", "$<$97"), ("all", "all")]


def simulated_tables(sim):
    """Table S3: (a) recall by read identity, origin and novelty band; (b) expected
    accuracy and reported organellar share for two sample compositions."""
    a_rows, b_rows = [], []
    band = {r["source_id"]: r["band"] for r in load(os.path.join(sim, "sources", "sources.tsv"))}
    for profile, label in PROFILES:
        path = os.path.join(sim, "bench", profile, "calls.tsv")
        if not os.path.exists(path):
            continue
        rows = load(path)
        src = {r["read_id"]: r["source_id"] for r in load(os.path.join(sim, "reads", f"{profile}.info.tsv"))}
        if a_rows:
            a_rows.append(r"\midrule")
            b_rows.append(r"\midrule")
        first = label
        rate = {}
        for o in ORIGINS:
            for b, b_label in BANDS:
                rs = [r for r in rows if r["truth"] == o and (b == "all" or band.get(src.get(r["read_id"])) == b)]
                if not rs:
                    continue
                vals = [f"{100 * sum(r[c] == o for r in rs) / len(rs):.1f}" for c, _ in METHODS]
                a_rows.append(f"{first} & {o if b == '99-100' else ''} & {b_label} & {num(len(rs))} & "
                              + " & ".join(vals) + r" \\")
                first = ""
            rs = [r for r in rows if r["truth"] == o]
            rate[o] = {c: Counter(r[c] for r in rs) for c, _ in METHODS}
            for c, _ in METHODS:
                rate[o][c] = {k: v / len(rs) for k, v in rate[o][c].items()}
        first = label
        for mix, w in MIXES:
            true_org = round(100 * (w["mitochondria"] + w["chloroplast"]))
            acc = [f"{100 * sum(w[o] * rate[o][c].get(o, 0) for o in ORIGINS):.1f}" for c, _ in METHODS]
            org = [f"{100 * sum(w[o] * (rate[o][c].get('mitochondria', 0) + rate[o][c].get('chloroplast', 0)) for o in ORIGINS):.1f}"
                   for c, _ in METHODS]
            b_rows.append(f"{first} & {mix}, accuracy & " + " & ".join(acc) + r" \\")
            b_rows.append(f" & {mix}, organellar (true {true_org}) & " + " & ".join(org) + r" \\")
            first = ""
    head = [r" & & \multicolumn{2}{c}{Metaxa2} & \multicolumn{3}{c}{SSUplex} \\"]
    out = [r"\begin{table}[htbp]", r"\centering", r"\scriptsize",
           r"\caption*{\textbf{(a) Recall (\%) by read identity, true origin and novelty band.} Band: "
           r"identity of the source to its closest Metaxa2 database sequence; ``all'' pools the bands. "
           r"Each source gave two reads.}",
           r"\begin{tabular}{lllrrrrrr}", r"\toprule",
           r" & & & & \multicolumn{2}{c}{Metaxa2} & \multicolumn{3}{c}{SSUplex} \\",
           r"\cmidrule(lr){5-6}\cmidrule(lr){7-9}",
           r"Read identity & Origin & Band (\%) & Reads & HMM & BLAST & sum & mean & count \\", r"\midrule"]
    out += a_rows + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    out += [r"\begin{table}[htbp]", r"\centering", r"\small",
            r"\caption*{\textbf{(b) Expected accuracy (\%) and reported organellar share (\%) for two "
            r"sample compositions.} Computed from the per-origin call rates in (a), which do not depend on "
            r"composition. Bacteria-dominated: 90\% bacteria, 2\% archaea, 2\% eukaryota, 3\% mitochondria, "
            r"3\% chloroplast (true organellar share 6\%). Plant-root-like: 55\% bacteria, 1\% archaea, "
            r"2\% eukaryota, 7\% mitochondria, 35\% chloroplast (42\%), close to Metaxa2's BLAST-assisted "
            r"calls on the rice-root sample.}",
            r"\begin{tabular}{llrrrrr}", r"\toprule"] + head + [
            r"\cmidrule(lr){3-4}\cmidrule(lr){5-7}",
            r"Read identity & Composition, measure & HMM & BLAST & sum & mean & count \\", r"\midrule"]
    out += b_rows + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accuracy", required=True)
    ap.add_argument("--sim")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    data = {k: load(os.path.join(a.accuracy, k, "calls.tsv")) for k, _ in LABELLED}
    rice = load(os.path.join(a.accuracy, "rice", "calls.tsv"))

    with open(os.path.join(a.out, "tableS1.tex"), "w") as f:
        for k, label in LABELLED:
            f.write(confusion(data[k], label))
        f.write(rice_table(rice))

    lines = [r"\begin{table}[htbp]", r"\centering", r"\scriptsize",
             r"\begin{tabular}{lllrlrrrrrr}", r"\toprule",
             r"\multicolumn{3}{c}{Origin chosen by} & & & \multicolumn{3}{c}{Origin chosen by mean} "
             r"& \multicolumn{3}{c}{Origin chosen by sum} \\",
             r"\cmidrule(lr){1-3}\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
             r"mean & sum & count & Reads & Truth (or Metaxa2 BLAST) & Regions & Mean & Sum "
             r"& Regions & Mean & Sum \\", r"\midrule"]
    for k, label in LABELLED:
        lines += disagreement_table(label, data[k], "truth") + [r"\midrule"]
    lines += disagreement_table("Rice-root sample (Metaxa2 BLAST call shown)", rice, "Metaxa2_final")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    with open(os.path.join(a.out, "tableS2.tex"), "w") as f:
        f.write("\n".join(lines))

    if a.sim:
        with open(os.path.join(a.out, "tableS3.tex"), "w") as f:
            f.write(simulated_tables(a.sim))
    print(f"wrote supplementary tables to {a.out}")


if __name__ == "__main__":
    main()
