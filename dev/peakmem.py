#!/usr/bin/env python3
"""Run a command and report its wall time, CPU time and peak memory.

Peak memory is reported two ways:
  tree_rss_mb  largest total resident set size of the command plus all of its
               descendants, sampled every --interval seconds (never reported
               below max_rss_mb, which very short runs can otherwise miss)
  max_rss_mb   largest single process in the tree, from getrusage(); this is
               the number /usr/bin/time reports, and it undercounts tools that
               run several processes at the same time

Both include the memory of the Python process that launches the command (the
kernel keeps the fork's high-water mark across exec), so values below about
20 MB measure this script rather than the tool. Benchmark runs are far larger.

Prints one tab-separated line: wall_s tree_rss_mb max_rss_mb cpu_s exit_code

Usage: peakmem.py [--interval 0.05] [--log FILE] -- command [args ...]
Requires psutil (pip install psutil, or conda install psutil).
"""
import argparse
import resource
import subprocess
import sys
import time

try:
    import psutil
except ImportError:
    sys.exit("peakmem.py needs psutil: pip install psutil (or conda install psutil)")


def tree_rss(proc):
    try:
        members = [proc] + proc.children(recursive=True)
    except psutil.NoSuchProcess:
        return 0
    total = 0
    for p in members:
        try:
            total += p.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interval", type=float, default=0.05)
    ap.add_argument("--log", help="file for the command's stdout and stderr")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        ap.error("no command given")

    log = open(a.log, "w") if a.log else subprocess.DEVNULL
    start = time.perf_counter()
    child = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT if a.log else subprocess.DEVNULL)
    proc = psutil.Process(child.pid)
    peak = 0
    while child.poll() is None:
        peak = max(peak, tree_rss(proc))
        time.sleep(a.interval)
    wall = time.perf_counter() - start

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    max_rss = usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024  # bytes
    cpu = usage.ru_utime + usage.ru_stime
    # The tree total can never be below its largest process; on very short runs
    # the sampler may miss the peak, so fall back to that lower bound.
    peak = max(peak, max_rss)
    print(f"{wall:.2f}\t{peak / 1e6:.1f}\t{max_rss / 1e6:.1f}\t{cpu:.1f}\t{child.returncode}")


if __name__ == "__main__":
    main()
