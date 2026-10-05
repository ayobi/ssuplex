//! Thin wrapper around `nhmmer` from HMMER 3.
//!
//! `nhmmer` (not `hmmsearch`) is the right tool for nucleotide profile vs
//! nucleotide read search: it scans *both* strands automatically, which matters
//! for real sequencing reads (e.g. ONT) that arrive in mixed orientation. The
//! wrapper:
//! 1. Discovers per-origin `.hmm` profile files in a directory
//! 2. Writes the in-memory reads to temporary FASTA files in chunks
//! 3. Invokes `nhmmer --tblout` per (chunk, profile)
//! 4. Parses the `--tblout` table into [`HmmHit`] structs

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{anyhow, bail, Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use tempfile::NamedTempFile;

use crate::cli::Args;
use crate::io::IndexedRead;
use crate::origin::Origin;

/// A single hit row from `nhmmer --tblout`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HmmHit {
    pub read_id: String,
    pub origin: Origin,
    /// Query model (conserved region) name, e.g. `B05`. Distinct regions are
    /// counted as independent evidence during origin resolution.
    pub region: String,
    pub score: f64,
    pub evalue: f64,
    /// Envelope coordinates on the query sequence (1-based, inclusive).
    pub env_from: usize,
    pub env_to: usize,
}

/// Origin profile loaded from disk.
#[derive(Debug, Clone)]
pub struct OriginProfile {
    pub origin: Origin,
    pub hmm_path: PathBuf,
}

/// Verify `nhmmer` is present and runnable.
///
/// Distinguishes not-on-`$PATH` from present-but-broken (e.g. an install missing
/// a shared library) and surfaces the captured stderr so the cause is visible.
pub fn ensure_nhmmer_available() -> Result<()> {
    let output = match Command::new("nhmmer").arg("-h").output() {
        Ok(o) => o,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            bail!("nhmmer not found on PATH — install HMMER 3 (e.g. `conda install -c bioconda hmmer`) and ensure nhmmer is on your PATH");
        }
        Err(e) => return Err(e).context("could not invoke nhmmer"),
    };

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stderr = stderr.trim();
        let code = output
            .status
            .code()
            .map(|c| c.to_string())
            .unwrap_or_else(|| "terminated by signal".to_string());
        bail!(
            "nhmmer is on PATH but failed to run (exit {code}){}",
            if stderr.is_empty() {
                String::new()
            } else {
                format!(": {stderr}")
            }
        );
    }
    Ok(())
}

/// Discover origin HMM profiles in `dir`.
///
/// Expects files named `{origin}.hmm` — e.g. `archaea.hmm`. Missing profiles
/// are silently skipped (you can run with just bacteria + archaea, for example).
pub fn load_origin_profiles(dir: &Path) -> Result<Vec<OriginProfile>> {
    if !dir.exists() {
        bail!("HMM directory does not exist: {}", dir.display());
    }

    let mut profiles = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) != Some("hmm") {
            continue;
        }
        if let Some(origin) = Origin::from_hmm_path(&path) {
            profiles.push(OriginProfile {
                origin,
                hmm_path: path,
            });
        }
    }

    if profiles.is_empty() {
        bail!(
            "no recognised origin .hmm files in {} (expected e.g. bacteria.hmm, archaea.hmm)",
            dir.display()
        );
    }

    profiles.sort_by_key(|p| p.origin.slug());
    Ok(profiles)
}

/// Scan all reads against all profiles, returning every qualifying hit.
///
/// Reads are chunked (`args.chunk_size`) and each (chunk × profile) pair is
/// dispatched in parallel via rayon. The progress bar tracks completed jobs.
pub fn scan_all(
    reads: &[IndexedRead],
    profiles: &[OriginProfile],
    args: &Args,
) -> Result<Vec<HmmHit>> {
    let chunks: Vec<&[IndexedRead]> = reads.chunks(args.chunk_size).collect();
    let jobs: Vec<(usize, &OriginProfile)> = chunks
        .iter()
        .enumerate()
        .flat_map(|(i, _)| profiles.iter().map(move |p| (i, p)))
        .collect();

    let pb = ProgressBar::new(jobs.len() as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
        )
        .unwrap()
        .progress_chars("##-"),
    );

    let results: Result<Vec<Vec<HmmHit>>> = jobs
        .par_iter()
        .map(|(chunk_idx, profile)| {
            let hits = scan_chunk(chunks[*chunk_idx], profile, args)?;
            pb.inc(1);
            Ok(hits)
        })
        .collect();

    pb.finish_and_clear();

    let mut flat: Vec<HmmHit> = results?.into_iter().flatten().collect();
    flat.shrink_to_fit();
    Ok(flat)
}

/// Run `nhmmer` for a single (chunk, profile) pair.
fn scan_chunk(chunk: &[IndexedRead], profile: &OriginProfile, args: &Args) -> Result<Vec<HmmHit>> {
    // Write chunk to a temporary FASTA.
    let mut query_file = NamedTempFile::new().context("failed to create temp query FASTA")?;
    for read in chunk {
        writeln!(query_file, ">{}", read.id)?;
        // Wrap sequence at 80 chars — HMMER is tolerant but it's the polite default.
        for line in read.sequence.as_bytes().chunks(80) {
            query_file.write_all(line)?;
            query_file.write_all(b"\n")?;
        }
    }
    query_file.flush()?;

    let tbl = NamedTempFile::new().context("failed to create temp tblout file")?;

    let output = Command::new("nhmmer")
        .arg("--tblout")
        .arg(tbl.path())
        .arg("-E")
        .arg(format!("{:e}", args.evalue))
        .arg("--noali")
        .arg("--cpu")
        .arg("1") // intra-process parallelism handled by rayon at the job level
        .arg(&profile.hmm_path)
        .arg(query_file.path())
        .stdout(Stdio::null()) // the human-readable report is large; we read --tblout instead
        .output()
        .with_context(|| format!("failed to invoke nhmmer for {}", profile.origin))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stderr = stderr.trim();
        let code = output
            .status
            .code()
            .map(|c| c.to_string())
            .unwrap_or_else(|| "signal".to_string());
        bail!(
            "nhmmer failed for origin {} (exit {}) on profile {}: {}",
            profile.origin,
            code,
            profile.hmm_path.display(),
            if stderr.is_empty() {
                "(no stderr captured)"
            } else {
                stderr
            }
        );
    }

    parse_tblout(tbl.path(), profile.origin)
}

/// Parse a `nhmmer --tblout` file into hits.
///
/// nhmmer tabular columns (whitespace-separated; the HMM is the *query*, the
/// reads are the *target*):
///   1: target name   (read_id)
///   3: query name     (conserved-region HMM model)
///   9: env from
///  10: env to         (env from > env to on the minus strand)
///  12: strand         (+/-)
///  13: E-value
///  14: score
///
/// nhmmer scans both strands, so a read in either orientation is found. On the
/// minus strand the envelope coordinates are reported high→low; we normalise to
/// (min, max) so downstream clipping is orientation-agnostic.
fn parse_tblout(path: &Path, origin: Origin) -> Result<Vec<HmmHit>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);

    // Best-scoring hit per (read, region) for this origin profile.
    let mut best: HashMap<(String, String), HmmHit> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 15 {
            continue;
        }

        let read_id = fields[0].to_string();
        let region = fields[2].to_string();
        let a: usize = fields[8]
            .parse()
            .map_err(|e| anyhow!("could not parse env from '{}': {e}", fields[8]))?;
        let b: usize = fields[9]
            .parse()
            .map_err(|e| anyhow!("could not parse env to '{}': {e}", fields[9]))?;
        let (env_from, env_to) = (a.min(b), a.max(b));
        let evalue: f64 = fields[12]
            .parse()
            .map_err(|e| anyhow!("could not parse E-value '{}': {e}", fields[12]))?;
        let score: f64 = fields[13]
            .parse()
            .map_err(|e| anyhow!("could not parse score '{}': {e}", fields[13]))?;

        let hit = HmmHit {
            read_id: read_id.clone(),
            origin,
            region: region.clone(),
            score,
            evalue,
            env_from,
            env_to,
        };

        best.entry((read_id, region))
            .and_modify(|existing| {
                if hit.score > existing.score {
                    *existing = hit.clone();
                }
            })
            .or_insert(hit);
    }

    Ok(best.into_values().collect())
}
