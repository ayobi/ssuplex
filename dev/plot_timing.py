#!/usr/bin/env python3
"""Summarise dev/run_timing.sh output and draw the performance figure.

Usage: plot_timing.py <timing.tsv> <out_prefix>
Writes <out_prefix>.summary.tsv (median, min and max per tool and read count)
and, if matplotlib is available, <out_prefix>.pdf and <out_prefix>.png. Only
measured points are plotted.
"""
import csv
import statistics
import sys
from collections import defaultdict

LABELS = {"ssuplex": "SSUplex", "metaxa2": "Metaxa2"}


def main():
    tsv, prefix = sys.argv[1], sys.argv[2]
    runs = defaultdict(list)
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["exit"] != "0":
                print(f"warning: {r['tool']} {r['reads']} rep {r['rep']} exited {r['exit']}; excluded")
                continue
            runs[(r["tool"], int(r["reads"]))].append(r)

    rows = []
    for (tool, n), rs in sorted(runs.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        row = {"tool": tool, "reads": n, "reps": len(rs)}
        for k in ("wall_s", "tree_rss_mb", "max_rss_mb", "cpu_s"):
            v = [float(r[k]) for r in rs]
            row[k] = statistics.median(v)
            row[k + "_min"], row[k + "_max"] = min(v), max(v)
        rows.append(row)
    cols = list(rows[0].keys())
    with open(prefix + ".summary.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(f"{row[c]:.2f}" if isinstance(row[c], float) else str(row[c]) for c in cols) + "\n")

    by = {(r["tool"], r["reads"]): r for r in rows}
    print(f"{'reads':>9}  {'SSUplex s':>10}  {'Metaxa2 s':>10}  {'speed-up':>8}  "
          f"{'SSUplex MB':>10}  {'Metaxa2 MB':>10}  {'mem ratio':>9}")
    for n in sorted({r["reads"] for r in rows}):
        s, m = by.get(("ssuplex", n)), by.get(("metaxa2", n))
        if s and m:
            print(f"{n:>9}  {s['wall_s']:>10.1f}  {m['wall_s']:>10.1f}  {m['wall_s'] / s['wall_s']:>7.2f}x  "
                  f"{s['tree_rss_mb']:>10.0f}  {m['tree_rss_mb']:>10.0f}  {s['tree_rss_mb'] / m['tree_rss_mb']:>9.2f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("matplotlib not available; wrote the summary table only")

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    for tool, marker in (("ssuplex", "o"), ("metaxa2", "s")):
        pts = [r for r in rows if r["tool"] == tool]
        if not pts:
            continue
        x = [r["reads"] for r in pts]
        for ax, key, scale in ((axes[0], "wall_s", 60.0), (axes[1], "tree_rss_mb", 1000.0)):
            y = [r[key] / scale for r in pts]
            lo = [(r[key] - r[key + "_min"]) / scale for r in pts]
            hi = [(r[key + "_max"] - r[key]) / scale for r in pts]
            ax.errorbar(x, y, yerr=[lo, hi], marker=marker, capsize=2, label=LABELS.get(tool, tool))
    axes[0].set_ylabel("Wall-clock time (min)")
    axes[1].set_ylabel("Peak memory, whole process tree (GB)")
    for ax, title in zip(axes, ("(a) Runtime", "(b) Peak memory")):
        ax.set_xscale("log")
        ax.set_xlabel("Reads")
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(prefix + ".pdf")
    fig.savefig(prefix + ".png", dpi=200)
    print(f"wrote {prefix}.summary.tsv, {prefix}.pdf, {prefix}.png")


if __name__ == "__main__":
    main()
