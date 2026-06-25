# Real-HMM benchmark: SSUplex vs Metaxa2

The smoke test (`dev/smoke.sh`) proves the pipeline mechanics on synthetic
motifs. This is the next step: run SSUplex with **Metaxa2's own SSU HMM
profiles** on a real read set, run Metaxa2 on the same reads, and compare
per-read origin calls. That tells us how faithfully the v0.1 reimplementation
tracks Metaxa2 — and exactly where it diverges.

## Prerequisites

```bash
conda install -c bioconda metaxa2 hmmer   # Metaxa2 + HMMER
cargo --version                            # modern stable toolchain
```

Find Metaxa2's database directory (it contains `SSU/HMMs/`). On a bioconda
install it is usually under `$CONDA_PREFIX/share/metaxa2*/metaxa2_db`, or
alongside the script:

```bash
which metaxa2
ls "$(dirname "$(readlink -f "$(which metaxa2)")")"      # look for metaxa2_db
find "$CONDA_PREFIX" -type d -name metaxa2_db 2>/dev/null
```

## A test dataset

Metaxa2 ships a small labelled test file — per its manual, ~50 SSU entries
(ten of each origin), some LSU, and a handful of non-SSU/non-LSU negatives.
That is the ideal first input: known origins, and built-in negatives to probe
false positives. Locate it inside the Metaxa2 install (look for `test`/`test_data`).
After that, scale up to a real eDNA set — a mock community, or one of the
acidophile/sediment runs where archaeal and chloroplast separation gets stressed.

## Run it

One command does assemble → Metaxa2 → SSUplex → compare:

```bash
dev/run_benchmark.sh \
  --db    /path/to/metaxa2_db \
  --reads /path/to/reads.fasta \
  --out   bench_out \
  --threads 8
```

Or step through manually:

```bash
# 1. assemble Metaxa2's profiles into SSUplex's layout
dev/assemble_metaxa2_hmms.sh /path/to/metaxa2_db bench_out/hmms
#    ^ INSPECT its output: it prints every .hmm it found and flags any it
#      could not bucket by origin. The filenames inside metaxa2_db vary by
#      version; if anything is listed as un-bucketed, sort it by hand before
#      trusting the run.

# 2. Metaxa2
metaxa2 -i reads.fasta -o bench_out/metaxa2/run --cpu 8 -g ssu

# 3. SSUplex with the same profiles
cargo build --release
target/release/ssuplex -i reads.fasta -o bench_out/ssuplex/sample \
  --hmm-dir bench_out/hmms -t 8

# 4. compare
dev/compare_origins.py \
  --ssuplex bench_out/ssuplex/sample.extraction.tsv \
  --metaxa2    bench_out/metaxa2/run.extraction.results \
  --reads      reads.fasta --out bench_out/concordance
```

## A note on Metaxa2's profiles

Metaxa2's SSU profiles use single-letter codes: A=archaea, B=bacteria,
C=chloroplast, E=eukaryota, M=mitochondria, and N=mitozoa (metazoan/animal
mitochondrial SSU). The assembler folds **N into mitochondria** so SSUplex's
five-origin scheme still detects animal mito reads; the comparison maps any
"mitozoa" label to mitochondria to match. The assembler prefers `SSU/HMMs`, so
pointing `--db` at the metaxa2_db root is fine — it won't sweep in the LSU or
COI profile sets.

## Reading the output

`compare_origins.py` prints a confusion matrix (rows = Metaxa2, cols =
SSUplex) and five buckets:

- **agree / disagree** — among reads both tools assigned to a concrete origin.
- **SSUplex-only** — SSUplex called an origin, Metaxa2 detected nothing.
- **Metaxa2-only** — Metaxa2 called an origin, SSUplex left it unclassified.
- **neither** — both declined.

### Expect divergence — it is informative, not a bug

100% concordance is *not* the target, because v0.1 deliberately simplifies
Metaxa2's classifier:

1. **The two-domain rule.** Metaxa2 requires at least two conserved SSU domains
   to hit before it classifies a read (per the Metaxa2 manual), which drives its
   false-positive rate near zero. SSUplex v0.1 takes the single best-scoring
   hit. So SSUplex will classify some reads Metaxa2 leaves undetected — that
   is the **SSUplex-only** bucket, and those are the calls to scrutinise as
   likely false positives. Closing this gap (a `--min-domains` criterion) is a
   natural v0.x feature; this benchmark quantifies how much it matters.
2. **Best-score-wins vs reliability scoring.** Metaxa2 weighs multiple regions
   and a reliability score; SSUplex uses a single best bit score across the
   per-origin profiles. Most disagreements will cluster on the genuinely hard
   pairs — bacterial vs chloroplast 16S, archaeal vs mitochondrial — which is
   exactly where to look.
3. **Coordinate agreement.** Beyond the origin label, compare envelope
   coordinates for reads both tools extracted (join on read_id via the CSV,
   then diff `env_from`/`env_to` against Metaxa2's reported positions).

### A reasonable first-pass success bar

- agreement among both-concrete reads is high (say >95%) on the labelled test set;
- disagreements concentrate on the known-hard origin pairs, not scattered randomly;
- the SSUplex-only bucket is small and explicable (sub-threshold/short reads).

Record per-run wall-clock and peak RSS too — `run_benchmark.sh` captures both
via `/usr/bin/time -v` into `*.time.log`, which feeds the throughput numbers in
`benches/README.md`.

## Length-stress benchmark (short reads)

Full-length reference SSU is the easy regime — both tools score near-perfect.
The regime that matters for real Illumina (V3–V4 ~250–450 bp) and partial ONT
reads is short reads that cover only one or two conserved regions, where
`--min-domains` and mean-score behave very differently. Simulate it from the
same labelled set by fragmenting the reference sequences (truth labels are
unchanged), isolating read length as the only variable:

```bash
# 300 bp windows instead of full ~1.5 kb
python3 dev/make_truth_set.py --db <metaxa2_db> --out /tmp/bench300 \
    --per-origin 20 --fragment-length 300

bash dev/run_benchmark.sh --db <metaxa2_db> \
    --reads /tmp/bench300.reads.fasta --truth /tmp/bench300.truth.tsv \
    --out bench300 --threads 8
```

What to look for as length drops: both tools' accuracy falls and the
unclassified fraction rises (fewer conserved regions are covered). `run_benchmark.sh`
forwards `--min-domains`, so the recall/precision trade is one command each:

```bash
bash dev/run_benchmark.sh --db <db> --reads /tmp/bench300.reads.fasta \
    --truth /tmp/bench300.truth.tsv --out bench300_d1 --min-domains 1
bash dev/run_benchmark.sh --db <db> --reads /tmp/bench300.reads.fasta \
    --truth /tmp/bench300.truth.tsv --out bench300_d2 --min-domains 2
```

Measured result at 300 bp: `--min-domains 2` lost ~11 correct calls to the
unclassified bucket and fixed zero errors, so **1 is the right default** for
short reads; raise it only when false positives dominate and reads span
multiple regions. Sweep lengths (150 / 300 / 600 / full) to map the curve.

### Isolating sequencing error

`make_truth_set.py --error-rate R` substitutes each base with probability R
(reference reads only), so you can separate noise from length:

```bash
python3 dev/make_truth_set.py --db <db> --out /tmp/bench300e2 \
    --per-origin 20 --fragment-length 300 --error-rate 0.02
```

Real Illumina is ~0.1–1% substitution; ONT is higher and indel-dominated (this
knob does substitutions only, so it under-models ONT — a floor, not a ceiling).

## Classifier design is locked at mean-primary

The origin ranking (mean per-region bit score, region count as tie-break) is
validated and intentional — see the DESIGN NOTE in `src/origin.rs`. Do not
switch to sum/count-weighted scoring without re-running this benchmark; `sum`
biases toward mitochondria (the largest profile after the N=mitozoa fold) and
regresses the result.

## Origin-ranking metric: evidence for the `sum` default

`--rank` selects how candidate origins are scored. The default (`sum`) was
chosen by running both benchmarks across all three metrics:

```text
                 reference (500, clean, 5 origins)   ZymoBIOMICS ONT (5000, real, bacteria)
  sum  (default)            95.6%                                95.9%
  mean                      96.8%                                43.0%
  count                     78.0%                                95.7%
```

`mean` wins marginally on clean full-length reads but collapses on real noisy
ONT (a wrong organelle origin beats the true origin on mean over fewer regions).
`count` over-collapses to bacteria on clean multi-origin data. `sum` is near-best
on both. The remaining `sum` errors on the reference set are the intrinsically
ambiguous bacteria-vs-mitochondria (alpha-proteobacterial SSU) reads, not a
scoring defect. Re-run BOTH benchmarks before changing the default.
