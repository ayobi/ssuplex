# dev/ — local test harness

Scripts for exercising the full pipeline end-to-end on synthetic fixtures.
These are **synthetic conserved motifs**, not real rRNA — they validate the
mechanics (HMM detection → origin call → envelope clip → per-origin output →
unclassified bucket), not biological accuracy.

## Requirements

- `cargo` (modern stable)
- `python3`
- HMMER on `PATH` (`hmmbuild` + `hmmsearch`) — `conda install -c bioconda hmmer`

## Run

```bash
# build release + run all checks
bash dev/smoke.sh

# or test an already-built binary
bash dev/smoke.sh target/debug/ssuplex
```

## What it checks

6 reads in (2 bacterial, 2 archaeal, 2 pure noise), then asserts:

- bacteria / archaea / unclassified FASTAs each have the expected record count
- the two no-hit reads are captured in `unclassified.fasta` (not silently dropped)
- summary counts are honest (`classified = 4`, `unclassified = 2`)
- no-hit reads are recorded in the extraction TSV
- envelope clipping lands on the predictable coordinates (e.g. `read_bact_1`
  has a 60 bp left flank → `coords=61-150`)

## Files

- `gen_fixtures.py <outdir>` — writes `bact.afa`, `arch.afa`, `reads.fasta`
  (deterministic, seed 42)
- `smoke.sh [binary]` — generates fixtures, builds the HMMs, runs the tool,
  asserts outputs
