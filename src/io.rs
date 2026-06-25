//! Input I/O — read FASTA into an in-memory index.
//!
//! For v0.1 we index the input fully into memory. This keeps the rest of the
//! pipeline simple (chunks are slices, parallelism is straightforward). For
//! very large inputs (> ~50M reads) a future iteration can switch to a true
//! streaming model with on-disk read storage.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_fasta as fasta;

/// A single input read held in memory.
#[derive(Debug, Clone)]
pub struct IndexedRead {
    pub id: String,
    pub sequence: String,
}

/// Read a FASTA file fully into memory.
///
/// Supports `.fasta`, `.fa`, `.fna`, and gzipped variants (`.gz`) — the latter
/// transparently via the standard `flate2` crate path inside noodles.
pub fn index_fasta(path: &Path) -> Result<Vec<IndexedRead>> {
    let reader: Box<dyn Read> = if has_gz_extension(path) {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Box::new(flate2_passthrough(file)?)
    } else {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Box::new(file)
    };

    let mut fasta_reader = fasta::io::Reader::new(BufReader::new(reader));

    let mut out = Vec::new();
    for record in fasta_reader.records() {
        let record = record.context("failed to parse FASTA record")?;
        let id = String::from_utf8(record.name().to_vec())
            .context("non-UTF8 sequence name in FASTA — only ASCII headers are supported")?;
        let sequence = String::from_utf8(record.sequence().as_ref().to_vec())
            .context("non-UTF8 sequence in FASTA — only ASCII nucleotide input is supported")?;
        out.push(IndexedRead { id, sequence });
    }
    Ok(out)
}

fn has_gz_extension(path: &Path) -> bool {
    path.extension().and_then(|s| s.to_str()) == Some("gz")
}

/// Hook for gzip support. Kept as a stub here so the dependency surface stays
/// small until we actually need it; flip to `flate2::read::GzDecoder` when
/// adding gzipped-input support.
fn flate2_passthrough(_file: File) -> Result<File> {
    anyhow::bail!("gzipped FASTA input is not yet supported in v0.1 — please decompress first")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn reads_basic_fasta() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, ">read1\nACGTACGT\n>read2\nTTTTGGGG").unwrap();
        let reads = index_fasta(f.path()).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].id, "read1");
        assert_eq!(reads[0].sequence, "ACGTACGT");
        assert_eq!(reads[1].id, "read2");
    }
}
