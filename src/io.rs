//! Input I/O — read FASTA into an in-memory index.
//!
//! The whole input is held in memory, which keeps the pipeline simple (chunks
//! are slices and parallelism is straightforward). Peak memory therefore
//! grows with the number and length of input reads.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use noodles_fasta as fasta;

/// A single input read held in memory.
#[derive(Debug, Clone)]
pub struct IndexedRead {
    pub id: String,
    pub sequence: String,
}

/// Read a FASTA file fully into memory.
///
/// Supports `.fasta`, `.fa`, `.fna`, and their gzipped (`.gz`) variants, which
/// are decompressed on the fly with `flate2`.
pub fn index_fasta(path: &Path) -> Result<Vec<IndexedRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn Read> = if has_gz_extension(path) {
        Box::new(GzDecoder::new(file))
    } else {
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

    #[test]
    fn reads_gzipped_fasta() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let f = tempfile::Builder::new()
            .suffix(".fasta.gz")
            .tempfile()
            .unwrap();
        let mut enc = GzEncoder::new(f.reopen().unwrap(), Compression::default());
        enc.write_all(b">g1\nACGTACGT\n>g2\nTTTTGGGG\n").unwrap();
        enc.finish().unwrap();
        let reads = index_fasta(f.path()).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].id, "g1");
        assert_eq!(reads[0].sequence, "ACGTACGT");
        assert_eq!(reads[1].sequence, "TTTTGGGG");
    }
}
