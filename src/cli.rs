use std::path::PathBuf;

use clap::{Parser, ValueEnum};

/// How to rank candidate origins from a read's per-region HMM hits.
///
/// For each origin, a read's hits (at most one per conserved-region profile)
/// are summarised as a region count, a summed bit score and a mean bit score.
/// These mirror Metaxa2's `--selection_priority` options `domains`, `sum` and
/// `score` (Metaxa2's default is `score`, the mean), except that Metaxa2
/// divides its sum by the number of profiles in the origin's set and SSUplex
/// does not. How the three compare on benchmark data is documented in
/// `dev/BENCHMARK.md`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum RankMetric {
    /// Mean per-region bit score; region count breaks ties
    Mean,
    /// Summed per-region bit score; region count breaks ties (default)
    Sum,
    /// Number of matched conserved regions; mean bit score breaks ties
    Count,
}

/// Fast preprocessor for rRNA-marker eDNA workflows. Extracts SSU rRNA
/// sequences from environmental reads and sorts them by origin (bacterial,
/// archaeal, eukaryotic, mitochondrial, chloroplast) for downstream
/// taxonomic classification. Reimplements the SSU extraction and
/// origin-assignment step of Metaxa2 in Rust.
#[derive(Parser, Debug, Clone)]
#[command(name = "ssuplex", version, author, about)]
pub struct Args {
    /// Input FASTA of sequencing reads or contigs (optionally gzip-compressed).
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    pub input: PathBuf,

    /// Output prefix. Files will be written as `{prefix}.{origin}.fasta`,
    /// `{prefix}.extraction.tsv`, and `{prefix}.summary.txt`.
    #[arg(short = 'o', long = "output", value_name = "PREFIX")]
    pub output: PathBuf,

    /// Directory containing origin-specific HMM profiles
    /// (`bacteria.hmm`, `archaea.hmm`, `eukaryota.hmm`,
    ///  `mitochondria.hmm`, `chloroplast.hmm`).
    ///
    /// Expects a directory of per-origin profiles built from the Metaxa2
    /// SSU database; see the README ("HMM profiles") for how to assemble it.
    #[arg(long = "hmm-dir", value_name = "DIR", default_value = "hmms")]
    pub hmm_dir: PathBuf,

    /// Number of worker threads [default: number of available CPUs]
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = num_cpus_default(),
        hide_default_value = true
    )]
    pub threads: usize,

    /// HMMER E-value cutoff for `nhmmer`.
    #[arg(short = 'E', long = "evalue", default_value_t = 1e-5)]
    pub evalue: f64,

    /// Minimum mean per-region bit score the chosen origin must reach;
    /// otherwise the read is reported as `unclassified`.
    #[arg(long = "min-score", default_value_t = 0.0)]
    pub min_score: f64,

    /// Minimum number of conserved-region profiles that must match the chosen
    /// origin; otherwise the read is reported as `unclassified`. (Metaxa2's
    /// default is 2, accepting single-region hits only at E-value 1e-10 or
    /// lower.)
    #[arg(long = "min-domains", default_value_t = 1)]
    pub min_domains: usize,

    /// Statistic used to rank candidate origins for each read: `sum` (summed
    /// per-region bit score, the default), `mean` (mean per-region bit score),
    /// or `count` (number of matched regions). See the README section
    /// "Origin ranking".
    #[arg(long = "rank", value_enum, default_value_t = RankMetric::Sum)]
    pub rank: RankMetric,

    /// Number of sequences per HMMER batch. Larger batches reduce subprocess
    /// overhead; smaller batches improve memory locality and progress
    /// reporting granularity.
    #[arg(long = "chunk-size", default_value_t = 10_000)]
    pub chunk_size: usize,

    /// Keep extracted reads as full input sequences (no clipping). By default,
    /// SSUplex clips each read to the SSU coordinate envelope from HMM hits.
    #[arg(long = "no-clip")]
    pub no_clip: bool,

    /// Also write `{prefix}.scores.tsv` (per read and matched origin:
    /// n_regions, mean, sum, best_evalue) and `{prefix}.regions.tsv` (every
    /// per-region hit with score and coordinates), for diagnostics.
    #[arg(long = "debug-scores")]
    pub debug_scores: bool,

    /// Replace output files left by a previous run with the same prefix.
    /// Without this flag SSUplex stops instead of overwriting existing results.
    #[arg(long = "force")]
    pub force: bool,

    /// Increase logging verbosity. Repeat (`-vv`) for trace output.
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count)]
    pub verbose: u8,
}

fn num_cpus_default() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_rank_is_sum() {
        let args = Args::try_parse_from(["ssuplex", "-i", "reads.fasta", "-o", "out"]).unwrap();
        assert_eq!(args.rank, RankMetric::Sum);
    }
}
