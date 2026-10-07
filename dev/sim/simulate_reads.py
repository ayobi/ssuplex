#!/usr/bin/env python3
"""Simulate full-length long reads from labelled sources with Badread's error models.

Each source sequence (and each negative) yields --reads-per-source reads per profile.
A read covers the whole source, as in amplicon sequencing, and comes from either
strand with equal probability. Errors are added with Badread's empirical models
(Wick 2019, https://github.com/rrwick/Badread), used as a library so that reads
stay full length:

  ont_r9     nanopore2020 model, identity beta distribution: mean 92%, sd 3%, max 98%
  ont_r10    nanopore2023 model, identity beta distribution: mean 98.5%, sd 1%, max 99.8%
  hifi_like  pacbio2021 model, qscores normal: mean Q30, sd 3 (about 99.9% identity)

Outputs per profile in --out:
  <profile>.fasta       reads (IDs <profile>_000001, ...)
  <profile>.truth.tsv   read_id, true_origin (or "(negative)"), for dev/score_origins.py
  <profile>.info.tsv    read_id, source_id, strand, target_identity, read_identity
"""
import argparse
import io
import multiprocessing as mp
import os
import random
import sys

import numpy as np

PROFILES = {
    "ont_r9": ("nanopore2020", 92.0, 3.0, 98.0),
    "ont_r10": ("nanopore2023", 98.5, 1.0, 99.8),
    "hifi_like": ("pacbio2021", 30.0, 3.0, None),
}
COMP = str.maketrans("ACGTN", "TGCAN")
_STATE = {}


def read_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_labels(path, column):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index(column)
        return {line.split("\t")[0]: line.rstrip("\n").split("\t")[idx] for line in f if line.strip()}


def init_worker(profile):
    from badread.error_model import ErrorModel
    from badread.identities import Identities
    from badread.qscore_model import QScoreModel
    model, mean, sd, mx = PROFILES[profile]
    quiet = io.StringIO()
    _STATE.update(error=ErrorModel(model, quiet), qscore=QScoreModel(model, quiet),
                  identities=Identities(mean, sd, mx, quiet))


def simulate_one(task):
    from badread.simulate import sequence_fragment
    import edlib
    seed, source_id, seq = task
    random.seed(seed)
    np.random.seed(seed % (2 ** 32))
    strand = "+" if random.random() < 0.5 else "-"
    template = seq if strand == "+" else seq.translate(COMP)[::-1]
    target = _STATE["identities"].get_identity()
    read, _, _, _ = sequence_fragment(template, target, _STATE["error"], _STATE["qscore"])
    dist = edlib.align(template, read)["editDistance"]
    return source_id, strand, target, 1 - dist / max(len(template), len(read)), read


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", required=True, help="directory written by build_sources.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--profiles", default=",".join(PROFILES))
    ap.add_argument("--reads-per-source", type=int, default=2)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    origin = read_labels(os.path.join(a.sources, "sources.tsv"), "origin")
    seqs = list(read_fasta(os.path.join(a.sources, "sources.fasta")))
    seqs += list(read_fasta(os.path.join(a.sources, "negatives.fasta")))
    for sid, _ in seqs:
        origin.setdefault(sid, "(negative)")

    for p_index, profile in enumerate(a.profiles.split(",")):
        if profile not in PROFILES:
            sys.exit(f"unknown profile {profile}; choose from {', '.join(PROFILES)}")
        tasks = [(a.seed * 1_000_003 + p_index * 100_000_007 + i * 101 + r, sid, s)
                 for i, (sid, s) in enumerate(seqs) for r in range(a.reads_per_source)]
        random.Random(f"{a.seed}-{profile}").shuffle(tasks)   # mix origins in the output
        print(f">> {profile}: simulating {len(tasks)} reads with {a.threads} processes", file=sys.stderr)
        with mp.Pool(a.threads, initializer=init_worker, initargs=(profile,)) as pool, \
             open(os.path.join(a.out, f"{profile}.fasta"), "w") as fa, \
             open(os.path.join(a.out, f"{profile}.truth.tsv"), "w") as truth, \
             open(os.path.join(a.out, f"{profile}.info.tsv"), "w") as info:
            truth.write("read_id\ttrue_origin\n")
            info.write("read_id\tsource_id\tstrand\ttarget_identity\tread_identity\n")
            for n, (sid, strand, target, ident, read) in enumerate(
                    pool.imap(simulate_one, tasks, chunksize=50), 1):
                rid = f"{profile}_{n:06d}"
                fa.write(f">{rid}\n{read}\n")
                truth.write(f"{rid}\t{origin[sid]}\n")
                info.write(f"{rid}\t{sid}\t{strand}\t{target:.4f}\t{ident:.4f}\n")
                if n % 20000 == 0:
                    print(f"   {n} reads", file=sys.stderr)
        print(f">> {profile}: wrote {len(tasks)} reads", file=sys.stderr)


if __name__ == "__main__":
    main()
