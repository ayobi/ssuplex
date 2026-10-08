#!/usr/bin/env python3
"""Summarise dev/run_timing.sh output and draw the performance figure.

Usage: plot_timing.py <timing.tsv> <out_prefix> [--full <timing_full.tsv>]
                      [--metaxa2-label LABEL]

Writes <out_prefix>.summary.tsv (median, minimum and maximum per tool and read
count) and, if matplotlib is available, <out_prefix>.pdf and <out_prefix>.png.
Only measured points are plotted, on log-log axes so that every input size is
readable and linear scaling appears as a straight line. Error bars span the
minimum and maximum over replicates.

--full adds Metaxa2's default pipeline (with BLAST classification) from a
second run_timing.sh output made with --m2-mode full; only its Metaxa2 rows
are used.
"""
import argparse
import csv
import statistics
import sys
from collections import defaultdict

KEYS = ("wall_s", "tree_rss_mb", "max_rss_mb", "cpu_s")


def summarise(path, tools=None):
    runs = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if tools and r["tool"] not in tools:
                continue
            if r["exit"] != "0":
                print(f"warning: {r['tool']} {r['reads']} rep {r['rep']} exited {r['exit']}; excluded")
                continue
            runs[(r["tool"], int(r["reads"]))].append(r)
    rows = []
    for (tool, n), rs in sorted(runs.items()):
        row = {"tool": tool, "reads": n, "reps": len(rs)}
        for k in KEYS:
            v = [float(r[k]) for r in rs]
            row[k] = statistics.median(v)
            row[k + "_min"], row[k + "_max"] = min(v), max(v)
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("timing")
    ap.add_argument("prefix")
    ap.add_argument("--full", help="run_timing.sh output made with --m2-mode full")
    ap.add_argument("--metaxa2-label", default="Metaxa2, extraction only (-x T)",
                    help="legend label for the Metaxa2 rows of <timing.tsv>")
    a = ap.parse_args()

    series = {"ssuplex": summarise(a.timing, {"ssuplex"}), "metaxa2": summarise(a.timing, {"metaxa2"})}
    if a.full:
        series["metaxa2_full"] = [dict(r, tool="metaxa2_full") for r in summarise(a.full, {"metaxa2"})]

    rows = [r for s in series.values() for r in s]
    cols = list(rows[0].keys())
    with open(a.prefix + ".summary.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(f"{row[c]:.2f}" if isinstance(row[c], float) else str(row[c]) for c in cols) + "\n")

    by = {(r["tool"], r["reads"]): r for r in rows}
    print(f"{'reads':>9}  {'SSUplex s':>10}  {'Metaxa2 s':>10}  {'speed-up':>8}  "
          f"{'SSUplex GB':>10}  {'Metaxa2 GB':>10}  {'mem ratio':>9}  {'full s':>8}  {'vs full':>8}")
    for n in sorted({r["reads"] for r in rows}):
        s, m, full = by.get(("ssuplex", n)), by.get(("metaxa2", n)), by.get(("metaxa2_full", n))
        if not s:
            continue
        line = f"{n:>9}  {s['wall_s']:>10.1f}"
        if m:
            line += (f"  {m['wall_s']:>10.1f}  {m['wall_s'] / s['wall_s']:>7.2f}x  {s['tree_rss_mb'] / 1000:>10.2f}"
                     f"  {m['tree_rss_mb'] / 1000:>10.2f}  {s['tree_rss_mb'] / m['tree_rss_mb']:>9.2f}")
        if full:
            line += f"  {full['wall_s']:>8.0f}  {full['wall_s'] / s['wall_s']:>7.1f}x"
        print(line)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("matplotlib not available; wrote the summary table only")

    style = {
        "ssuplex": ("SSUplex", "o", "#0072B2"),
        "metaxa2": (a.metaxa2_label, "s", "#D55E00"),
        "metaxa2_full": ("Metaxa2, full pipeline with BLAST", "^", "#CC79A7"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    for key, pts in series.items():
        if not pts:
            continue
        label, marker, colour = style[key]
        x = [r["reads"] for r in pts]
        for ax, k, scale in ((axes[0], "wall_s", 60.0), (axes[1], "tree_rss_mb", 1000.0)):
            y = [r[k] / scale for r in pts]
            lo = [(r[k] - r[k + "_min"]) / scale for r in pts]
            hi = [(r[k + "_max"] - r[k]) / scale for r in pts]
            # open markers for the full pipeline, whose memory coincides with -x T
            face = "white" if key == "metaxa2_full" else colour
            ax.errorbar(x, y, yerr=[lo, hi], marker=marker, color=colour, markerfacecolor=face,
                        capsize=2, label=label)
    axes[0].set_ylabel("Wall-clock time (min)")
    axes[1].set_ylabel("Peak memory, whole process tree (GB)")
    for ax, title in zip(axes, ("(a) Run time", "(b) Peak memory")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Reads")
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(True, which="both", alpha=0.3)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels), fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(a.prefix + ".pdf")
    fig.savefig(a.prefix + ".png", dpi=200)
    print(f"wrote {a.prefix}.summary.tsv, {a.prefix}.pdf, {a.prefix}.png")


if __name__ == "__main__":
    main()
