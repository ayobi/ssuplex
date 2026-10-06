# dev/

Scripts for testing and benchmarking SSUplex. See [BENCHMARK.md](BENCHMARK.md)
for the full protocol.

| Script | Purpose |
| ------ | ------- |
| `assemble_metaxa2_hmms.sh` | Build the five per-origin profile files from a Metaxa2 database |
| `fetch_ena_reads.sh` | Download the first N reads of an ENA/SRA run as FASTA |
| `make_truth_set.py` | Sample a labelled reference set from the Metaxa2 SSU database |
| `run_real_data.sh` | Run the accuracy benchmark on all labelled and real datasets |
| `run_benchmark.sh` | Run Metaxa2 and SSUplex (mean, sum, count) on one dataset and score them |
| `score_origins.py` | Score origin calls against truth and against Metaxa2 |
| `compare_rank_rules.py` | Compare origin-ranking rules offline on the hits from a benchmark run |
| `run_timing.sh` | Measure run time and peak memory across input sizes |
| `peakmem.py` | Run one command and report wall time and process-tree peak memory |
| `plot_timing.py` | Summarise timing results and draw the performance figure |
| `smoke.sh`, `gen_fixtures.py` | End-to-end test on synthetic profiles and reads |
