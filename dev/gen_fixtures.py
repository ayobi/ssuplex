#!/usr/bin/env python3
"""Generate deterministic test fixtures for the SSUplex smoke test.

Writes, into <outdir>:
  bact.afa     aligned FASTA training set for a synthetic "bacterial" motif
  arch.afa     aligned FASTA training set for a synthetic "archaeal" motif
  reads.fasta  6 reads: 2 bacterial, 2 archaeal, 2 pure-noise

These are *synthetic* conserved motifs, not real 16S — they exercise the
pipeline mechanics (HMM detection -> best-score origin call -> envelope clip
-> per-origin output -> unclassified bucket), not biological accuracy.

Usage:  python3 dev/gen_fixtures.py <outdir>
"""
import os
import random
import sys

SEED = 42
BASES = "ACGT"


def rnd(n, rng):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(s, k, rng):
    s = list(s)
    for _ in range(k):
        i = rng.randrange(len(s))
        s[i] = rng.choice(BASES)
    return "".join(s)


def write_afa(path, seqs, tag):
    with open(path, "w") as f:
        for i, s in enumerate(seqs):
            f.write(f">{tag}_{i}\n{s}\n")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: gen_fixtures.py <outdir>")
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(SEED)

    # Two clearly distinct conserved 90-bp motifs.
    bact_motif = rnd(90, rng)
    arch_motif = rnd(90, rng)
    arch_motif = "".join(
        c if rng.random() > 0.7 else rng.choice(BASES) for c in arch_motif
    )

    # Training alignments: 8 sequences each, lightly mutated, equal length
    # (no indels, so they are already "aligned" for hmmbuild).
    def aln(motif, n, k):
        return [mutate(motif, k, rng) for _ in range(n)]

    write_afa(os.path.join(outdir, "bact.afa"), aln(bact_motif, 8, 4), "bact")
    write_afa(os.path.join(outdir, "arch.afa"), aln(arch_motif, 8, 4), "arch")

    # Reads: embed a (mutated) motif in random flanks of known length, so the
    # reported envelope coordinates are predictable.
    def read(name, motif, lflank, rflank, k=3):
        return f">{name}\n{rnd(lflank, rng)}{mutate(motif, k, rng)}{rnd(rflank, rng)}\n"

    with open(os.path.join(outdir, "reads.fasta"), "w") as f:
        f.write(read("read_bact_1", bact_motif, 60, 50))  # motif at 61..150
        f.write(read("read_bact_2", bact_motif, 20, 120))  # motif at 21..110
        f.write(read("read_arch_1", arch_motif, 100, 30))  # motif at 101..190
        f.write(read("read_arch_2", arch_motif, 5, 5, k=2))  # motif at 6..95
        f.write(f">read_junk_1\n{rnd(200, rng)}\n")  # pure noise -> no hit
        f.write(f">read_junk_2\n{rnd(150, rng)}\n")  # pure noise -> no hit

    print(f"wrote fixtures to {outdir}/ (bact.afa, arch.afa, reads.fasta)")


if __name__ == "__main__":
    main()
