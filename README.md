# SSUplex

**Fast extraction and origin-sorting of small-subunit rRNA (16S/18S) from environmental sequencing reads.**

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.07.02.736232-b31b1b)](https://doi.org/10.64898/2026.07.02.736232)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/ayobi/ssuplex/actions/workflows/ci.yml/badge.svg)](https://github.com/ayobi/ssuplex/actions/workflows/ci.yml)
[![Rust](https://img.shields.io/badge/rust-1.85%2B-b7410e?logo=rust&logoColor=white)](https://www.rust-lang.org/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://anaconda.org/bioconda/ssuplex)
[![GitHub release](https://img.shields.io/github/v/release/ayobi/ssuplex)](https://github.com/ayobi/ssuplex/releases)

SSUplex finds small-subunit rRNA in sequencing reads, assigns each read to one
of five origins (bacteria, archaea, eukaryota, mitochondria, chloroplast),
clips it to the SSU region and writes one FASTA file per origin, so that each
set can be passed to an appropriate taxonomic classifier. It reimplements the
HMM-based SSU extraction and origin-assignment step of
[Metaxa2](https://microbiology.se/software/metaxa2/) in Rust, using Metaxa2's
own SSU profiles. It does not perform taxonomic classification.

## Why sort reads by origin

Broadly conserved rRNA primers amplify more than the intended target: a 16S run
from plant or host material can contain many chloroplast and mitochondrial
reads, and 18S or cross-domain reads also appear. Sent directly to a bacterial
16S classifier, chloroplast reads tend to be assigned to *Cyanobacteria*,
mitochondrial reads to *Rickettsiales*, and eukaryotic reads to low-confidence
hits. SSUplex separates these reads first. For each read it:

1. **detects** SSU rRNA with profile HMMs of the gene's conserved regions,
   searching both strands with `nhmmer`, so reads in either orientation are
   found without a separate orientation step;
2. **assigns an origin** from the per-region scores of each origin's profiles;
3. **extracts** the SSU part of the read, clipped to the envelope of its hits.

## Installation

### Bioconda (recommended)

```bash
conda create -n ssuplex -c conda-forge -c bioconda ssuplex
conda activate ssuplex
ssuplex --version
```

This also installs HMMER, which provides `nhmmer`, the only runtime dependency.

### From source

You need Rust 1.85 or newer (install with [rustup](https://rustup.rs)) and
HMMER 3.1 or newer with `nhmmer` on your `PATH` (for example
`conda install -c bioconda hmmer` or `brew install hmmer`).

```bash
git clone https://github.com/ayobi/ssuplex.git
cd ssuplex
cargo install --locked --path .
ssuplex --version
```

`cargo install` builds a release binary and copies it to `~/.cargo/bin`, which
rustup adds to your `PATH`. `--locked` builds with exactly the dependency
versions recorded in `Cargo.lock`. If your shell reports
`ssuplex: command not found`, either add `~/.cargo/bin` to `PATH` or build in
place with `cargo build --release --locked` and run `./target/release/ssuplex`.

## HMM profiles

SSUplex needs one profile file per origin (`bacteria.hmm`, `archaea.hmm`,
`eukaryota.hmm`, `mitochondria.hmm`, `chloroplast.hmm`) in the directory given
by `--hmm-dir` (default `./hmms`). They are built from the Metaxa2 SSU
database, which is distributed with Metaxa2 under its own license, so they are
not bundled. A working Metaxa2 installation is not required, only its download:

```bash
curl -L -O https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz
tar xzf Metaxa2_2.2.3.tar.gz
dev/assemble_metaxa2_hmms.sh Metaxa2_2.2.3/metaxa2_db hmms
```

The script concatenates Metaxa2's per-origin profile files (A archaea,
B bacteria, C chloroplast, E eukaryota, and M plus N for mitochondria, where N
is Metaxa2's metazoan mitochondrial set) and prints how many profiles each file
contains. The equivalent commands are:

```bash
SRC=Metaxa2_2.2.3/metaxa2_db/SSU/HMMs
mkdir -p hmms
cat "$SRC/B.hmm"              > hmms/bacteria.hmm
cat "$SRC/A.hmm"              > hmms/archaea.hmm
cat "$SRC/E.hmm"              > hmms/eukaryota.hmm
cat "$SRC/M.hmm" "$SRC/N.hmm" > hmms/mitochondria.hmm
cat "$SRC/C.hmm"              > hmms/chloroplast.hmm
```

## Quick start

The Metaxa2 download includes `test.fasta`, 100 sequences with known content:
50 SSU sequences (ten from each origin) and 50 LSU sequences, which SSUplex
should leave unclassified. From the directory where you assembled `hmms/`:

```bash
ssuplex -i Metaxa2_2.2.3/test.fasta -o quickstart/test --hmm-dir hmms --threads 4
cat quickstart/test.summary.txt
```

The log should end with `classified 50 of 100 read(s) to an origin`. With
HMMER 3.4 and default settings, the 50 SSU sequences are assigned 12 to
bacteria, 10 each to archaea and eukaryota, and 9 each to mitochondria and
chloroplast: one mitochondrial and one chloroplast sequence end up as bacteria
(see [Origin ranking](#origin-ranking)). Running
again with the same `-o` prefix stops with an error instead of overwriting the
results; add `--force` to replace them.

## Output files

For `-o results/sample01`:

| File                          | Contents                                                     |
| ----------------------------- | ------------------------------------------------------------ |
| `sample01.bacteria.fasta`     | Reads assigned to bacteria, clipped to the SSU region       |
| `sample01.archaea.fasta`      | Reads assigned to archaea                                    |
| `sample01.eukaryota.fasta`    | Reads assigned to eukaryota (nuclear 18S)                    |
| `sample01.mitochondria.fasta` | Reads assigned to mitochondria                               |
| `sample01.chloroplast.fasta`  | Reads assigned to chloroplast                                |
| `sample01.unclassified.fasta` | Reads without a qualifying SSU hit, unclipped                |
| `sample01.extraction.tsv`     | Per read: origin, score, region count, E-value, coordinates  |
| `sample01.summary.txt`        | Run settings and counts per origin                           |
| `sample01.scores.tsv`         | With `--debug-scores`: per read and origin, region count, mean and summed bit score |
| `sample01.regions.tsv`        | With `--debug-scores`: every per-region hit                  |

A per-origin FASTA file is written only when at least one read is assigned to
that origin.

## Options

```text
$ ssuplex -h
Fast Rust preprocessor for rRNA-marker eDNA workflows: extracts SSU rRNA from environmental sequencing reads and sorts them by origin (bacterial, archaeal, eukaryotic, mitochondrial, chloroplast).

Usage: ssuplex [OPTIONS] --input <FASTA> --output <PREFIX>

Options:
  -i, --input <FASTA>              Input FASTA of sequencing reads or contigs (optionally gzip-compressed)
  -o, --output <PREFIX>            Output prefix. Files will be written as `{prefix}.{origin}.fasta`, `{prefix}.extraction.tsv`, and `{prefix}.summary.txt`
      --hmm-dir <DIR>              Directory containing origin-specific HMM profiles (`bacteria.hmm`, `archaea.hmm`, `eukaryota.hmm`, `mitochondria.hmm`, `chloroplast.hmm`) [default: hmms]
  -t, --threads <THREADS>          Number of worker threads [default: number of available CPUs]
  -E, --evalue <EVALUE>            HMMER E-value cutoff for `nhmmer` [default: 0.00001]
      --min-score <MIN_SCORE>      Minimum mean per-region bit score the chosen origin must reach; otherwise the read is reported as `unclassified` [default: 0]
      --min-domains <MIN_DOMAINS>  Minimum number of conserved-region profiles that must match the chosen origin; otherwise the read is reported as `unclassified`. (Metaxa2's default is 2, accepting single-region hits only at E-value 1e-10 or lower.) [default: 1]
      --rank <RANK>                Statistic used to rank candidate origins for each read: `sum` (summed per-region bit score, the default), `mean` (mean per-region bit score), or `count` (number of matched regions). See the README section "Origin ranking" [default: sum] [possible values: mean, sum, count]
      --chunk-size <CHUNK_SIZE>    Number of sequences per HMMER batch. Larger batches reduce subprocess overhead; smaller batches improve memory locality and progress reporting granularity [default: 10000]
      --no-clip                    Keep extracted reads as full input sequences (no clipping). By default, SSUplex clips each read to the SSU coordinate envelope from HMM hits
      --debug-scores               Also write `{prefix}.scores.tsv` (per read and matched origin: n_regions, mean, sum, best_evalue) and `{prefix}.regions.tsv` (every per-region hit with score and coordinates), for diagnostics
      --force                      Replace output files left by a previous run with the same prefix. Without this flag SSUplex stops instead of overwriting existing results
  -v, --verbose...                 Increase logging verbosity. Repeat (`-vv`) for trace output
  -h, --help                       Print help (see more with '--help')
  -V, --version                    Print version
```

## Origin ranking

For each read and origin, SSUplex keeps at most one hit per conserved-region
profile and summarises the hits as a region count, a summed bit score and a mean
bit score. `--rank` chooses which of these decides the origin:

| `--rank`        | Origin chosen by                  | Ties broken by    |
| --------------- | --------------------------------- | ----------------- |
| `sum` (default) | highest summed bit score          | number of regions |
| `mean`          | highest mean per-region bit score | number of regions |
| `count`         | most matched regions              | mean bit score    |

These mirror Metaxa2's `--selection_priority` options `sum`, `score` (Metaxa2's
default) and `domains`, with one difference: Metaxa2 divides its sum by the
number of profiles in the origin's set, and SSUplex does not. Because
`mitochondria.hmm` combines Metaxa2's M and N sets, a mitochondrial origin can
collect hits from both. A read is reported as unclassified when the chosen
origin matches fewer than `--min-domains` regions or its mean score is below
`--min-score`.

**Why `sum` is the default.** On noisy long reads, the true origin of a read
typically matches more conserved regions than a related origin does, but each
match scores a little lower. The mean ignores how many regions matched: on the
ZymoBIOMICS Nanopore mock community it assigned more than half of the bacterial
reads to chloroplast or mitochondria, and so does Metaxa2's own HMM-based step,
which also ranks by the mean. The sum takes breadth into account and was the
only statistic that stayed accurate on both clean reference sequences and noisy
Nanopore reads ([`dev/BENCHMARK.md`](dev/BENCHMARK.md)).

Its weakness is the opposite case: organellar reads that also match many
bacterial profiles can be assigned to bacteria. In the benchmarks this affected
about one in five full-length chloroplast reference sequences and more than a
third of the reads that Metaxa2 (with BLAST) called mitochondrial in a
plant-root sample. Metaxa2 resolves such reads with its BLAST step; with
SSUplex, expect some organellar reads in the bacterial output. Classifiers whose
reference databases include organelle sequences, such as SILVA, will label them
as chloroplast or mitochondria.

## Relation to Metaxa2

SSUplex uses Metaxa2's SSU profiles and reproduces its extraction and
HMM-based origin assignment. Metaxa2's default pipeline then runs a BLAST
classification against its reference database, and its final origin calls
also take those BLAST matches into account; `metaxa2 -x T` runs the extraction
step alone. SSUplex has no classification step and leaves taxonomy to a
dedicated classifier. Other differences: SSUplex searches with `nhmmer`
(Metaxa2 uses `hmmsearch` on both strands), ranks origins by summed rather
than mean score by default, uses an E-value cutoff of 1e-5 and one matching
region by default (Metaxa2 requires two domains, or one domain at E-value 1e-10
or lower), and scores Metaxa2's M and N mitochondrial sets together.

## Scope and limitations

- SSUplex is built as a preprocessing step for 16S/18S amplicon metabarcoding.
  It will run on shotgun reads, but that use has not been benchmarked, and
  because the whole input is held in memory, peak memory grows with input size.
- Organellar reads that also match bacterial profiles are the least certain
  calls; see [Origin ranking](#origin-ranking).
- The profiles come from the Metaxa2 database (built from SILVA release 111 and
  Mitozoa release 10), so lineages that are poorly represented there may be
  detected less reliably.

## Reproducing the benchmarks

[`dev/BENCHMARK.md`](dev/BENCHMARK.md) describes the datasets, the exact
commands, and how accuracy, run time and memory are measured. `dev/smoke.sh`
runs an end-to-end test on synthetic profiles and reads.

## Citation

O'Brien A, Vargas Botia J, Acuña I, Restovic F, Martínez P, Parada P. SSUplex:
fast, both-strand extraction and origin-sorting of small-subunit rRNA for
environmental DNA metabarcoding. bioRxiv (2026).
doi: [10.64898/2026.07.02.736232](https://doi.org/10.64898/2026.07.02.736232)

Please also cite Metaxa2, whose profiles and extraction approach SSUplex uses:

> Bengtsson-Palme J, Hartmann M, Eriksson KM, Pal C, Thorell K, Larsson DGJ,
> Nilsson RH (2015). Metaxa2: Improved identification and taxonomic
> classification of small and large subunit rRNA in metagenomic data.
> *Molecular Ecology Resources* 15(6):1403-1414.
> doi: [10.1111/1755-0998.12399](https://doi.org/10.1111/1755-0998.12399)

## License

MIT; see [LICENSE](LICENSE). The HMM profiles you build are derived from the
Metaxa2 database and remain under Metaxa2's license.
