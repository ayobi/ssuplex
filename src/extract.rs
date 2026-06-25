//! SSU subsequence extraction.
//!
//! Given per-read origin calls with envelope coordinates from `nhmmer`,
//! clip each read down to the SSU envelope. With `--no-clip`, the whole input
//! sequence is preserved (useful when downstream classifiers want to see
//! flanking context).

use std::collections::HashMap;

use anyhow::Result;

use crate::io::IndexedRead;
use crate::origin::{Origin, OriginCall};

/// A read after origin-aware extraction.
#[derive(Debug, Clone)]
pub struct ExtractedRead {
    pub id: String,
    pub sequence: String,
    /// Slice coordinates on the original read (1-based, inclusive).
    pub start: usize,
    pub end: usize,
}

/// Build extracted reads for every call. Reads classified as `Unclassified`
/// are not extracted (returned in a separate bucket by [`crate::summary`]).
pub fn extract_all(
    reads: &[IndexedRead],
    calls: &[OriginCall],
) -> Result<HashMap<String, ExtractedRead>> {
    let by_id: HashMap<&str, &IndexedRead> = reads.iter().map(|r| (r.id.as_str(), r)).collect();

    let mut out = HashMap::with_capacity(calls.len());
    for call in calls {
        // Unclassified reads have no SSU envelope to clip to; the summary
        // writer emits them whole into `unclassified.fasta` instead.
        if call.origin == Origin::Unclassified {
            continue;
        }
        let Some(read) = by_id.get(call.read_id.as_str()) else {
            continue;
        };

        let len = read.sequence.len();
        // Clamp coordinates into [1, len] and convert to 0-based half-open.
        let start_1 = call.env_from.max(1).min(len);
        let end_1 = call.env_to.max(start_1).min(len);
        let start_0 = start_1 - 1;
        let end_0 = end_1; // half-open

        let sequence = read.sequence[start_0..end_0].to_string();

        out.insert(
            call.read_id.clone(),
            ExtractedRead {
                id: call.read_id.clone(),
                sequence,
                start: start_1,
                end: end_1,
            },
        );
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::origin::Origin;

    fn read(id: &str, seq: &str) -> IndexedRead {
        IndexedRead {
            id: id.to_string(),
            sequence: seq.to_string(),
        }
    }

    fn call(id: &str, from: usize, to: usize) -> OriginCall {
        OriginCall {
            read_id: id.to_string(),
            origin: Origin::Bacteria,
            score: 100.0,
            n_regions: 1,
            evalue: 1e-20,
            env_from: from,
            env_to: to,
        }
    }

    #[test]
    fn clips_to_envelope() {
        let reads = vec![read("r1", "AAAACCCCGGGGTTTT")]; // len 16
        let calls = vec![call("r1", 5, 12)]; // CCCCGGGG
        let out = extract_all(&reads, &calls).unwrap();
        let r = &out["r1"];
        assert_eq!(r.sequence, "CCCCGGGG");
        assert_eq!(r.start, 5);
        assert_eq!(r.end, 12);
    }

    #[test]
    fn clamps_out_of_bounds() {
        let reads = vec![read("r1", "ACGT")];
        let calls = vec![call("r1", 0, 999)];
        let out = extract_all(&reads, &calls).unwrap();
        assert_eq!(out["r1"].sequence, "ACGT");
    }
}
