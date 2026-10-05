//! Output writers — per-origin FASTAs, extraction TSV, run summary.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::cli::Args;
use crate::extract::ExtractedRead;
use crate::hmm::HmmHit;
use crate::io::IndexedRead;
use crate::origin::{debug_score_table, Origin, OriginCall};
use crate::VERSION;

/// Write all output files for the run.
pub fn write_all(
    args: &Args,
    reads: &[IndexedRead],
    calls: &[OriginCall],
    extracted: &HashMap<String, ExtractedRead>,
) -> Result<()> {
    write_per_origin_fastas(args, reads, calls, extracted)?;
    write_extraction_tsv(args, calls)?;
    write_summary(args, reads, calls)?;
    Ok(())
}

/// Write the `--debug-scores` per-(read, origin) table to `{prefix}.scores.tsv`.
pub fn write_debug_scores(args: &Args, reads: &[IndexedRead], hits: &[HmmHit]) -> Result<()> {
    let mut p = args.output.clone();
    let stem = p
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "ssuplex".to_string());
    p.set_file_name(format!("{stem}.scores.tsv"));

    let f = File::create(&p).with_context(|| format!("failed to create {}", p.display()))?;
    let mut w = BufWriter::new(f);
    writeln!(w, "read_id\torigin\tn_regions\tmean\tsum\tbest_evalue")?;
    for r in debug_score_table(reads, hits) {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2e}",
            r.read_id, r.origin, r.n_regions, r.mean, r.sum, r.best_evalue
        )?;
    }
    Ok(())
}

/// Write the raw per-region HMM hits to `{prefix}.regions.tsv`
/// (read_id, origin, region, score, evalue, env_from, env_to), sorted by
/// read then origin then descending score. Lets us see exactly which conserved
/// regions each origin matched on a read — the input to the per-origin means.
pub fn write_debug_regions(args: &Args, hits: &[HmmHit]) -> Result<()> {
    let mut p = args.output.clone();
    let stem = p
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "ssuplex".to_string());
    p.set_file_name(format!("{stem}.regions.tsv"));

    let mut sorted: Vec<&HmmHit> = hits.iter().collect();
    sorted.sort_by(|a, b| {
        a.read_id
            .cmp(&b.read_id)
            .then(a.origin.slug().cmp(b.origin.slug()))
            .then(
                b.score
                    .partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });

    let f = File::create(&p).with_context(|| format!("failed to create {}", p.display()))?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "read_id\torigin\tregion\tscore\tevalue\tenv_from\tenv_to"
    )?;
    for h in sorted {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.2}\t{:.2e}\t{}\t{}",
            h.read_id, h.origin, h.region, h.score, h.evalue, h.env_from, h.env_to
        )?;
    }
    Ok(())
}

/// Every file a run can write for `args.output`, in a fixed order. Used to
/// refuse overwriting a previous run and, with `--force`, to clear its files
/// (including per-origin FASTAs and debug tables the new run may not write).
pub fn output_paths(args: &Args) -> Vec<PathBuf> {
    let mut paths: Vec<PathBuf> = Origin::ALL
        .iter()
        .copied()
        .chain(std::iter::once(Origin::Unclassified))
        .map(|o| fasta_path_for(&args.output, o))
        .collect();
    for suffix in ["extraction.tsv", "summary.txt", "scores.tsv", "regions.tsv"] {
        let mut p = args.output.clone();
        let stem = p
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "ssuplex".to_string());
        p.set_file_name(format!("{stem}.{suffix}"));
        paths.push(p);
    }
    paths
}

fn fasta_path_for(prefix: &Path, origin: Origin) -> PathBuf {
    let mut p = prefix.to_path_buf();
    let stem = p
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "ssuplex".to_string());
    p.set_file_name(format!("{stem}.{}.fasta", origin.slug()));
    p
}

fn write_per_origin_fastas(
    args: &Args,
    reads: &[IndexedRead],
    calls: &[OriginCall],
    extracted: &HashMap<String, ExtractedRead>,
) -> Result<()> {
    let reads_by_id: HashMap<&str, &IndexedRead> =
        reads.iter().map(|r| (r.id.as_str(), r)).collect();

    // One writer per origin, opened lazily.
    let mut writers: HashMap<Origin, BufWriter<File>> = HashMap::new();

    for call in calls {
        let path = fasta_path_for(&args.output, call.origin);
        let writer = match writers.get_mut(&call.origin) {
            Some(w) => w,
            None => {
                if let Some(parent) = path.parent() {
                    if !parent.as_os_str().is_empty() {
                        std::fs::create_dir_all(parent).ok();
                    }
                }
                let f =
                    File::create(&path).with_context(|| format!("creating {}", path.display()))?;
                writers.insert(call.origin, BufWriter::new(f));
                writers.get_mut(&call.origin).unwrap()
            }
        };

        let (seq, coords) = if args.no_clip || call.origin == Origin::Unclassified {
            // Emit the full input sequence: either the user opted out of
            // clipping, or the read is unclassified and has no SSU envelope.
            let read = reads_by_id
                .get(call.read_id.as_str())
                .ok_or_else(|| anyhow::anyhow!("read {} not found in index", call.read_id))?;
            (read.sequence.as_str().to_string(), (1, read.sequence.len()))
        } else {
            let ex = extracted
                .get(&call.read_id)
                .ok_or_else(|| anyhow::anyhow!("no extracted record for {}", call.read_id))?;
            (ex.sequence.clone(), (ex.start, ex.end))
        };

        writeln!(
            writer,
            ">{} origin={} score={:.1} coords={}-{}",
            call.read_id, call.origin, call.score, coords.0, coords.1
        )?;
        for line in seq.as_bytes().chunks(80) {
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
        }
    }

    for (_, mut w) in writers {
        w.flush()?;
    }
    Ok(())
}

fn write_extraction_tsv(args: &Args, calls: &[OriginCall]) -> Result<()> {
    let mut path = args.output.clone();
    let stem = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "ssuplex".to_string());
    path.set_file_name(format!("{stem}.extraction.tsv"));

    let f = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);

    writeln!(
        w,
        "read_id\torigin\tscore\tn_regions\tevalue\tenv_from\tenv_to"
    )?;
    for c in calls {
        writeln!(
            w,
            "{}\t{}\t{:.2}\t{}\t{:.2e}\t{}\t{}",
            c.read_id, c.origin, c.score, c.n_regions, c.evalue, c.env_from, c.env_to
        )?;
    }
    w.flush()?;
    Ok(())
}

fn write_summary(args: &Args, reads: &[IndexedRead], calls: &[OriginCall]) -> Result<()> {
    let mut path = args.output.clone();
    let stem = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "ssuplex".to_string());
    path.set_file_name(format!("{stem}.summary.txt"));

    let f = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);

    let total = reads.len();
    let mut by_origin: HashMap<Origin, usize> = HashMap::new();
    for c in calls {
        *by_origin.entry(c.origin).or_insert(0) += 1;
    }
    let unclassified = by_origin.get(&Origin::Unclassified).copied().unwrap_or(0);
    let classified = total.saturating_sub(unclassified);

    writeln!(w, "SSUplex v{} run summary", VERSION)?;
    writeln!(w, "===========================")?;
    writeln!(w, "input:          {}", args.input.display())?;
    writeln!(w, "output prefix:  {}", args.output.display())?;
    writeln!(w, "threads:        {}", args.threads)?;
    writeln!(w, "E-value cutoff: {:e}", args.evalue)?;
    writeln!(w, "min bit score:  {}", args.min_score)?;
    writeln!(w, "min domains:    {}", args.min_domains)?;
    writeln!(w)?;
    writeln!(w, "total input reads:  {}", total)?;
    writeln!(w, "classified reads:   {}", classified)?;
    writeln!(w, "unclassified reads: {}", unclassified)?;
    writeln!(w)?;
    writeln!(w, "Per-origin counts:")?;
    for origin in Origin::ALL
        .iter()
        .chain(std::iter::once(&Origin::Unclassified))
    {
        let n = by_origin.get(origin).copied().unwrap_or(0);
        writeln!(w, "  {:<14} {:>10}", origin.to_string(), n)?;
    }
    w.flush()?;
    Ok(())
}
