# SSUplex

**A fast, lightweight Rust preprocessor for rRNA-marker eDNA workflows.**

SSUplex extracts small-subunit ribosomal RNA (SSU; 16S / 18S) sequences from
environmental sequencing reads and sorts them by origin..bacterial, archaeal,
eukaryotic, mitochondrial, or chloroplast! soo each pile can be routed to the
appropriate downstream taxonomic classifier. It is a Rust reimplementation of
the SSU extraction and origin-classification logic of
[Metaxa2](https://microbiology.se/software/metaxa2/), built for modern
long-read (Oxford Nanopore, PacBio HiFi) and short-read amplicon data.

---

## What it does

A single PCR with "universal" rRNA primers pulls in sequences from every domain
of life, plus organelle markers that ride along on bacterial primers. Fed
straight into a bacterial 16S classifier, chloroplast reads get misassigned as
*Cyanobacteria*, mitochondrial reads as *Rickettsiales*, and 18S reads return
noise. SSUplex does the triage. For each read it:

1. **Detects**:: SSU rRNA using HMM profiles of the conserved subregions of the
   gene, scanning both strands via `nhmmer` (so reverse-complement reads..
   roughly half of any long-read run are detected without a separate
   orientation step!!).
2. **Classifies origin**:: as bacterial, archaeal, eukaryotic, mitochondrial, or
   chloroplast.
3. **Extracts**:: the SSU portion, clipping to the HMM coordinate envelope.

It writes one FASTA per origin, plus a per-read coordinate/score table and a run
summary. Downstream taxonomic classification is left to the appropriate tool
(EMU, QIIME 2, mothur, nf-core/ampliseq) operating on each origin-sorted FASTA
against the proper reference database.

## Why

Metaxa2 has been the standard tool for this triage step for over a decade, and
it remains accurate and capable. Buttt it is built on Perl orchestration around
HMMER and BLAST, was last released in 2021 (v2.2.3), and folds taxonomic
classification (a BLAST stage) into every run. SSUplex is a focused
reimplementation of just the extraction-and-origin step in Rust:

- **A single static binary**:: no Perl runtime, trivial to install and embed in
  pipelines or containers.
- **Faster at scale**:: comparable to Metaxa2's extraction on small inputs, and
  increasingly faster as read counts grow, reaching ~3.4× at 200,000 reads
  while using less memory and fewer CPU cores (see [Performance](#performance)).
- **Extraction-only by design**:: it stops at origin-sorted reads and hands off
  to whichever dedicated classifier your study uses, so you don't pay for a
  BLAST classification stage you may not want.

## Workflow

```
eDNA reads → SSUplex ──┬─ bacteria.fasta     → bacterial classifier   (SILVA / GTDB)
                       ├─ archaea.fasta      → archaeal classifier    (KSGP+ / SILVA)
                       ├─ eukaryota.fasta    → 18S classifier         (PR2 / SILVA-Euk)
                       ├─ chloroplast.fasta  → plant/algal classifier (PhytoREF)
                       ├─ mitochondria.fasta → metazoan classifier    (MIDORI / MitoFish)
                       └─ unclassified.fasta → diagnostics
```

## Installation

### From source

```bash
git clone https://github.com/ayobi/ssuplex.git
cd ssuplex
cargo build --release
./target/release/ssuplex --help
```

### Runtime dependency

SSUplex shells out to `nhmmer` from [HMMER 3](http://hmmer.org/) and expects it
on `$PATH`:

```bash
conda install -c bioconda hmmer
```

### HMM profiles

SSUplex needs one HMM profile per origin — `bacteria.hmm`, `archaea.hmm`,
`eukaryota.hmm`, `mitochondria.hmm`, `chloroplast.hmm` — in a directory passed
with `--hmm-dir` (default: `./hmms`). These are **not** bundled with SSUplex:
they are built from the [Metaxa2](https://microbiology.se/software/metaxa2/) SSU
profile database, which carries its own license, so you assemble them once from a
Metaxa2 download. You do not need a working Metaxa2 install — only its profiles.

Download and unpack the Metaxa2 package:

```bash
curl -L -O https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz
tar xzf Metaxa2_2.2.3.tar.gz
```

Assemble the five origin files into `hmms/` with the helper script. It reads
`<metaxa2_db>/SSU/HMMs/`, buckets the profiles by origin, folds Metaxa2's
metazoan-mitochondrial (`N`) set into `mitochondria`, and prints what it did:

```bash
dev/assemble_metaxa2_hmms.sh Metaxa2_2.2.3/metaxa2_db hmms
```

The script uses `mapfile` and associative arrays, so it needs **bash 4+**. macOS
ships bash 3.2 — install a newer one (`brew install bash`) and run it explicitly:
`$(brew --prefix)/bin/bash dev/assemble_metaxa2_hmms.sh Metaxa2_2.2.3/metaxa2_db hmms`.

If you prefer not to run the script, the five files are just concatenations of
Metaxa2's single-letter SSU profiles (`A`=archaea, `B`=bacteria, `C`=chloroplast,
`E`=eukaryota, `M`=mitochondria, `N`=metazoan mitochondria):

```bash
SRC=Metaxa2_2.2.3/metaxa2_db/SSU/HMMs
mkdir -p hmms
cat "$SRC/B.hmm"              > hmms/bacteria.hmm
cat "$SRC/A.hmm"              > hmms/archaea.hmm
cat "$SRC/E.hmm"              > hmms/eukaryota.hmm
cat "$SRC/M.hmm" "$SRC/N.hmm" > hmms/mitochondria.hmm
cat "$SRC/C.hmm"              > hmms/chloroplast.hmm
```

Verify you have five non-empty files before running SSUplex:

```bash
ls -la hmms/
```

## Quick start

```bash
ssuplex -i reads.fasta -o results/sample01 --threads 8
```

| File                          | Contents                                               |
| ----------------------------- | ------------------------------------------------------ |
| `sample01.bacteria.fasta`     | Reads classified as bacterial SSU (16S)                |
| `sample01.archaea.fasta`      | Archaeal SSU (16S)                                      |
| `sample01.eukaryota.fasta`    | Eukaryotic SSU (18S)                                    |
| `sample01.mitochondria.fasta` | Mitochondrial SSU                                       |
| `sample01.chloroplast.fasta`  | Chloroplast SSU (16S)                                   |
| `sample01.unclassified.fasta` | Reads with no qualifying SSU hit                        |
| `sample01.extraction.tsv`     | Per-read origin call, region count, score, coordinates |
| `sample01.summary.txt`        | Run statistics                                          |

`--rank {mean,sum,count}` selects the origin-ranking statistic (see
[Validation](#validation)). `--min-domains N` sets the minimum conserved regions
for a call (default 1). `--debug-scores` additionally writes per-read,
per-origin score and per-region hit tables for diagnostics.

## Performance

Both tools were run extraction-only (SSUplex default; `metaxa2 -x T`), both
scanning both strands, both with 12 threads (`--cpu 12 --multi_thread T`), on
full-length ONT 16S reads (ZymoBIOMICS mock, SRR10391201) streamed from ENA, on
a single workstation. Times are wall-clock; the 5,000-read point was run in
triplicate with <2% variance.

| Reads   | SSUplex  | Metaxa2 `-x T` | Speedup | SSUplex peak RAM | Metaxa2 peak RAM |
| ------- | -------- | -------------- | ------- | ---------------- | ---------------- |
| 5,000   | 7.6 s    | 8.9 s          | 1.2×    | 56 MB            | 72 MB            |
| 25,000  | 18.7 s   | 45.5 s         | 2.4×    | 298 MB           | 328 MB           |
| 100,000 | 58.9 s   | 3 m 03 s       | 3.1×    | 946 MB           | 1.25 GB          |
| 200,000 | 1 m 49 s | 6 m 07 s       | 3.4×    | 1.58 GB          | 2.45 GB          |

At small inputs the two are comparable; the gap widens with scale, reaching
~3.4× at 200,000 reads (and still widening), while SSUplex uses ~35% less memory, 
and does so occupying fewer CPU cores (~4 vs ~6.5). This 50k–500k read range
is typical of per-sample amplicon sequencing. Note this compares extraction
only... Metaxa2's full pipeline additionally runs BLAST classification, which
SSUplex omits by design.

**Threading and memory notes.** SSUplex parallelises across the five origin
profiles, so multithreading delivers roughly a 3–4× speedup and saturates around
four threads, requesting more does not help! Memory scales with read count
(sequences are indexed in memory), which is comfortable at amplicon scale but
would be the limiting factor for very large (multi-million-read) metagenomes;
streamed input is on the roadmap.

## Validation

SSUplex is benchmarked against Metaxa2 with a reproducible harness (`dev/`) on
three public datasets: a full-length SSU reference set sampled from the Metaxa2
database, the ZymoBIOMICS ONT mock (SRR10391201; all-bacterial truth), and a
real plant-root ONT 16S sample (PRJNA992961, rice endophytes; scored as
concordance with Metaxa2, since it has no per-read truth).

The origin-ranking statistic (`--rank`) involves a genuine, measured trade-off:

| `--rank`         | Reference (clean) | Zymo ONT (bacteria) | Rice ONT (multi-origin)            |
| ---------------- | ----------------- | ------------------- | ---------------------------------- |
| `mean` (default) | 96.8%             | 43.0%               | tracks Metaxa2 on mito/chloroplast |
| `sum`            | 95.6%             | 95.9%               | misroutes plant mitochondria       |
| `count`          | 78.0%             | 95.7%               | over-collapses to bacteria         |

**Guidance:** the default `mean` reproduces Metaxa2 most faithfully on
full-length/clean reads and on organelle-rich samples. For bacteria-dominated
noisy long-read samples (mock communities, many gut/soil ONT datasets) where
organelles are absent, `--rank sum` is markedly more accurate. Choose by sample
type.

## Scope and limits

SSUplex is strong at:

- SSU detection and extraction, on both strands, from long and short reads.
- Distinguishing chloroplast, archaeal, and eukaryotic SSU, reliably across
  clean and noisy data.
- Bacterial classification on clean/full-length reads, and on noisy
  bacteria-dominated samples with `--rank sum`.

The **bacteria-vs-mitochondria** boundary is its lower-confidence edge :(
Mitochondrial SSU is alpha-proteobacterial in origin, so on noisy reads the two
are genuinely close, and the winner depends on region-inclusion details that no
single HMM-score statistic resolves cleanly. Metaxa2 disambiguates these reads
with a second BLAST classification stage; SSUplex does not (by design.. that
stage is the main runtime cost). Treat mitochondrial vs bacterial calls on noisy
data as lower-confidence, and prefer the default `mean` when mitochondrial
signal matters.

## Reproducing the benchmarks

See `dev/BENCHMARK.md`. In brief: `dev/make_truth_set.py` builds labelled read
sets from the Metaxa2 DB (with optional fragmentation and error injection);
`dev/run_benchmark.sh` runs Metaxa2 and SSUplex and prints a confusion matrix
via `dev/compare_origins.py`. Public ONT data is streamed from ENA (e.g.
`SRR10391201` for Zymo, `SRR25243163` for rice).

## Roadmap

- [x] v0.1: parallel SSU extraction, both-strand (`nhmmer`), origin classification, benchmark harness
- [ ] updated archaeal HMM panel (GTDB / KSGP+) for DPANN / Asgardarchaeota / Thermoplasmatota recovery
- [ ] ambiguous-read classification cross-check (lightweight mapping tie-breaker for the bacteria/mito boundary)
- [ ] streamed input for multi-million-read metagenomes
- [ ] long-read error model for ONT / PacBio HiFi
- [ ] LSU (23S / 28S) extraction
- [ ] V-region subextraction for hybrid long-read → short-read-classifier workflows

## Citation

Manuscript in preparation; cite the GitHub release for now. When publishing,
please also cite the original Metaxa2 paper whose extraction approach this tool
reimplements!

> Bengtsson-Palme J, Hartmann M, Eriksson KM, Pal C, Thorell K, Larsson DGJ, Nilsson RH (2015).
> Metaxa2: Improved identification and taxonomic classification of small and large subunit
> rRNA in metagenomic data. *Molecular Ecology Resources*, 15(6), 1403–1414.
> doi: [10.1111/1755-0998.12399](https://doi.org/10.1111/1755-0998.12399)

## License

MIT. See [LICENSE](LICENSE).
