use std::path::PathBuf;

use clap::{Parser, ValueEnum};

/// How to rank candidate origins from their per-region HMM hits.
///
/// Chosen empirically against two benchmarks (`dev/run_benchmark.sh`): a
/// full-length Metaxa2 reference set and the ZymoBIOMICS ONT mock, plus a real
/// plant-root ONT sample (rice, PRJNA992961) scored as concordance with Metaxa2:
///
/// ```text
///            reference    Zymo (bacteria)   rice (multi-origin, vs Metaxa2)
///   Mean       96.8%           43.0%        tracks Metaxa2 best (mito/chloro)
///   Sum        95.6%           95.9%        misroutes plant mitochondria
///   Count      78.0%           95.7%        over-collapses to bacteria
/// ```
///
/// `Mean` is the default: it tracks Metaxa2 most faithfully on multi-origin
/// data, where it recovers the plant-mitochondrial and chloroplast fractions
/// that `Sum`/`Count` push to bacteria. On bacteria-only samples (mocks, many
/// gut/soil datasets) `Sum` scores higher because piling on bacteria's many
/// conserved regions can't be wrong when everything is bacterial — so `--rank
/// sum` is offered for that case.
///
/// HONEST LIMIT: the bacteria-vs-mitochondria call cannot be made robustly from
/// HMM region scores alone (mito SSU is alpha-proteobacterial-derived). On noisy
/// reads the two are genuinely close and the winner depends on region-inclusion
/// details no single statistic resolves cleanly; Metaxa2 leans on a second BLAST
/// stage there, which SSUplex does not have. SSUplex is therefore strong
/// and fast on bacterial and chloroplast origin assignment, and flags the
/// bacteria/mitochondria boundary as its lower-confidence edge.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum RankMetric {
    /// Mean per-region bit score (region count breaks ties). Default; tracks
    /// Metaxa2 best on multi-origin data.
    Mean,
    /// Sum of per-region bit scores — total weight of evidence. Stronger on
    /// bacteria-dominated samples; misroutes organelle reads to bacteria.
    Sum,
    /// Number of distinct conserved regions matched (mean score breaks ties).
    /// Over-favours the broadest (bacterial) profile.
    Count,
}

/// Fast preprocessor for rRNA-marker eDNA workflows. Extracts SSU rRNA
/// sequences from environmental reads and sorts them by origin (bacterial,
/// archaeal, eukaryotic, mitochondrial, chloroplast) for downstream
/// taxonomic classification. A Rust reimplementation of Metaxa2.
#[derive(Parser, Debug, Clone)]
#[command(name = "ssuplex", version, author, about)]
pub struct Args {
    /// Input FASTA of sequencing reads or contigs.
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

    /// Number of worker threads.
    #[arg(short = 't', long = "threads", default_value_t = num_cpus_default())]
    pub threads: usize,

    /// HMMER E-value cutoff for `nhmmer`.
    #[arg(short = 'E', long = "evalue", default_value_t = 1e-5)]
    pub evalue: f64,

    /// Minimum bit score required for a read to be assigned to an origin.
    /// Reads with no hit above this score are reported as `unclassified`.
    #[arg(long = "min-score", default_value_t = 0.0)]
    pub min_score: f64,

    /// Minimum number of distinct conserved regions (HMM query models) that
    /// must match an origin for a read to be assigned to it. Metaxa2 uses 2,
    /// which keeps the false-positive rate low; the default here is 1 so that
    /// single-model profiles still work. Origins are ranked by mean per-region
    /// bit score, so raising this mainly suppresses weak single-region calls.
    #[arg(long = "min-domains", default_value_t = 1)]
    pub min_domains: usize,

    /// Statistic used to rank candidate origins for each read: `mean` (mean
    /// per-region bit score, the default), `sum` (total per-region bit score),
    /// or `count` (number of matched regions). `mean` tracks Metaxa2 best on
    /// multi-origin data; `sum` is stronger on bacteria-dominated samples. See
    /// `RankMetric` for the benchmark evidence.
    #[arg(long = "rank", value_enum, default_value_t = RankMetric::Mean)]
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

    /// Write a per-read, per-origin score table to `{prefix}.scores.tsv`
    /// (read_id, origin, n_regions, mean, sum, best_evalue) for every origin
    /// that matched at least one region. Diagnostic for tuning origin scoring.
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
