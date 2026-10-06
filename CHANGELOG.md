# Changelog

## Unreleased

### Changed
- The default origin-ranking statistic is now `sum` (previously `mean`). On
  noisy Nanopore reads the mean assigned many bacterial reads to chloroplast or
  mitochondria; see the README section "Origin ranking" and `dev/BENCHMARK.md`.
  Use `--rank mean` for the previous behaviour.
- SSUplex no longer overwrites the outputs of a previous run with the same
  prefix; it stops with an error. Use `--force` to replace them.
- Minimum supported Rust version is 1.85.

### Fixed
- Runs no longer abort when `nhmmer` cannot guess the alphabet of a chunk of
  very short or low-complexity reads; the alphabet is now passed explicitly.
- Help texts and documentation corrected: input is held in memory, `--min-score`
  applies to the chosen origin's mean score, gzip input is accepted, and
  `--debug-scores` also writes `regions.tsv`.

### Added
- `--force`.
- Benchmark harness in `dev/`: accuracy against truth and agreement with both
  Metaxa2 stages, run time and process-tree memory across input sizes, offline
  comparison of ranking rules, and scripts to fetch the public datasets.
- CI runs the README installation and quick start verbatim.

### Removed
- Unused `thiserror` dependency.

## 0.1.0 (2026-07-17)

Initial release.
