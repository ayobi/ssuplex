//! Origin classification.
//!
//! Given a set of HMM hits across several origin profiles (bacteria, archaea,
//! eukaryota, mitochondria, chloroplast), this module decides which origin a
//! read belongs to. Hits are aggregated per origin and ranked by **mean
//! per-region bit score**, with the number of matching conserved regions as a
//! tie-breaker, and an optional minimum-region floor (`min_domains`). Reads
//! with no qualifying origin are reported as [`Origin::Unclassified`].

use std::collections::HashMap;
use std::fmt;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::cli::RankMetric;
use crate::hmm::HmmHit;
use crate::io::IndexedRead;

/// Origin of an SSU rRNA read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Origin {
    Bacteria,
    Archaea,
    Eukaryota,
    Mitochondria,
    Chloroplast,
    Unclassified,
}

impl Origin {
    /// All five classifiable origins (excludes `Unclassified`).
    pub const ALL: [Origin; 5] = [
        Origin::Bacteria,
        Origin::Archaea,
        Origin::Eukaryota,
        Origin::Mitochondria,
        Origin::Chloroplast,
    ];

    /// Lowercase identifier used for file naming and HMM lookup.
    pub fn slug(&self) -> &'static str {
        match self {
            Origin::Bacteria => "bacteria",
            Origin::Archaea => "archaea",
            Origin::Eukaryota => "eukaryota",
            Origin::Mitochondria => "mitochondria",
            Origin::Chloroplast => "chloroplast",
            Origin::Unclassified => "unclassified",
        }
    }

    /// Parse from the slug used in HMM filenames (`bacteria.hmm`, etc.).
    pub fn from_slug(s: &str) -> Option<Origin> {
        match s.to_ascii_lowercase().as_str() {
            "bacteria" => Some(Origin::Bacteria),
            "archaea" => Some(Origin::Archaea),
            "eukaryota" | "eukaryote" => Some(Origin::Eukaryota),
            "mitochondria" | "mitochondrial" => Some(Origin::Mitochondria),
            "chloroplast" | "plastid" => Some(Origin::Chloroplast),
            _ => None,
        }
    }

    /// Infer origin from an HMM file path by its stem.
    pub fn from_hmm_path(path: &Path) -> Option<Origin> {
        path.file_stem()
            .and_then(|s| s.to_str())
            .and_then(Self::from_slug)
    }
}

impl fmt::Display for Origin {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.slug())
    }
}

/// Best-origin call for a single read.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OriginCall {
    pub read_id: String,
    pub origin: Origin,
    /// Mean per-region bit score across the conserved regions that matched the
    /// chosen origin (the metric origins are ranked by).
    pub score: f64,
    /// Number of distinct conserved regions that matched the chosen origin.
    pub n_regions: usize,
    pub evalue: f64,
    /// Envelope start on the read (1-based, inclusive). Spans all matched
    /// regions of the chosen origin (min over their envelope starts).
    pub env_from: usize,
    /// Envelope end on the read (1-based, inclusive). Spans all matched
    /// regions of the chosen origin (max over their envelope ends).
    pub env_to: usize,
}

/// Per-origin aggregate of a read's conserved-region hits.
#[derive(Clone)]
struct Agg {
    n: usize,
    sum_score: f64,
    best_evalue: f64,
    env_from: usize,
    env_to: usize,
}

impl Agg {
    fn mean(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.sum_score / self.n as f64
        }
    }
}

/// Aggregate hits per read, per origin: distinct conserved regions, summed
/// score, best E-value, and envelope span. Shared by [`resolve`] and
/// [`debug_score_table`] so they can never disagree on the aggregates.
fn aggregate(hits: &[HmmHit]) -> HashMap<&str, HashMap<Origin, Agg>> {
    let mut by_read: HashMap<&str, HashMap<Origin, Agg>> = HashMap::new();
    for hit in hits {
        let agg = by_read
            .entry(hit.read_id.as_str())
            .or_default()
            .entry(hit.origin)
            .or_insert(Agg {
                n: 0,
                sum_score: 0.0,
                best_evalue: f64::INFINITY,
                env_from: usize::MAX,
                env_to: 0,
            });
        agg.n += 1;
        agg.sum_score += hit.score;
        if hit.evalue < agg.best_evalue {
            agg.best_evalue = hit.evalue;
        }
        agg.env_from = agg.env_from.min(hit.env_from);
        agg.env_to = agg.env_to.max(hit.env_to);
    }
    by_read
}

/// One per-(read, origin) row for the `--debug-scores` diagnostic dump.
pub struct DebugScoreRow {
    pub read_id: String,
    pub origin: Origin,
    pub n_regions: usize,
    pub mean: f64,
    pub sum: f64,
    pub best_evalue: f64,
}

/// Emit one row per (read, origin) that matched at least one region, in input
/// read order (origins within a read sorted by slug for determinism). Lets the
/// per-origin scores SSUplex actually computes be analysed directly, rather
/// than inferred from Metaxa2's tallies.
pub fn debug_score_table(reads: &[IndexedRead], hits: &[HmmHit]) -> Vec<DebugScoreRow> {
    let by_read = aggregate(hits);
    let mut rows = Vec::new();
    for read in reads {
        let Some(origins) = by_read.get(read.id.as_str()) else {
            continue;
        };
        let mut ordered: Vec<(&Origin, &Agg)> = origins.iter().collect();
        ordered.sort_by_key(|(o, _)| o.slug());
        for (origin, agg) in ordered {
            rows.push(DebugScoreRow {
                read_id: read.id.clone(),
                origin: *origin,
                n_regions: agg.n,
                mean: agg.mean(),
                sum: agg.sum_score,
                best_evalue: agg.best_evalue,
            });
        }
    }
    rows
}

/// Resolve per-read origin calls.
///
/// Each read's hits are aggregated per origin: the number of distinct
/// conserved regions that matched, their summed/mean bit score, the best
/// E-value, and the envelope span. Origins are then ranked by the chosen
/// [`RankMetric`] (with a secondary key and origin slug as final, deterministic
/// tie-breaks). This replaces "single best hit wins": because chloroplast SSU is
/// bacterial-derived, one organelle region can out-score the best bacterial
/// region on a true bacterial read, even though bacteria wins on breadth.
///
/// DESIGN NOTE — choice of metric is empirical, not obvious, so it is exposed as
/// `--rank` and evaluated against `dev/run_benchmark.sh` on BOTH the full-length
/// reference set and the ZymoBIOMICS ONT mock:
///   * `Mean` scored 50/50 on clean full-length reference reads, but on noisy
///     real ONT reads it misroutes ~45% of bacteria to chloroplast/mito,
///     because a wrong organelle origin often has a higher mean on FEWER
///     regions while the true origin has MORE regions at a slightly lower mean.
///   * `Sum` and `Count` weight breadth and should recover those reads, but may
///     over-favour bacteria on true organelle reads.
///
/// Do not hard-code one metric or change the default without re-running both
/// benchmarks; record the numbers when you do.
///
/// The winning origin is accepted only if it matched at least `min_domains`
/// regions and its mean score is at least `min_score`; otherwise the read is
/// [`Origin::Unclassified`] (the best origin's aggregate is still recorded for
/// diagnostics). Reads with no hit at all become `Unclassified` with a sentinel.
/// Exactly one [`OriginCall`] is emitted per input read, in input order.
pub fn resolve(
    reads: &[IndexedRead],
    hits: &[HmmHit],
    min_score: f64,
    min_domains: usize,
    rank: RankMetric,
) -> Vec<OriginCall> {
    // read_id -> origin -> aggregate (shared with debug_score_table)
    let by_read = aggregate(hits);

    reads
        .iter()
        .map(|read| {
            let Some(origins) = by_read.get(read.id.as_str()) else {
                return OriginCall {
                    read_id: read.id.clone(),
                    origin: Origin::Unclassified,
                    score: 0.0,
                    n_regions: 0,
                    evalue: f64::INFINITY,
                    env_from: 0,
                    env_to: 0,
                };
            };

            // Rank origins by the chosen metric, with a metric-specific
            // secondary key and origin slug as the final deterministic
            // tie-break. (No score is NaN — values come from HMMER.)
            // Primary/secondary key per metric:
            //   Mean  -> (mean,      n   )
            //   Sum   -> (sum_score, n   )
            //   Count -> (n as f64,  mean)
            let key = |a: &Agg| -> (f64, f64) {
                match rank {
                    RankMetric::Mean => (a.mean(), a.n as f64),
                    RankMetric::Sum => (a.sum_score, a.n as f64),
                    RankMetric::Count => (a.n as f64, a.mean()),
                }
            };
            let (origin, agg) = origins
                .iter()
                .max_by(|(oa, a), (ob, b)| {
                    let (pa, sa) = key(a);
                    let (pb, sb) = key(b);
                    pa.partial_cmp(&pb)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(sa.partial_cmp(&sb).unwrap_or(std::cmp::Ordering::Equal))
                        .then(ob.slug().cmp(oa.slug()))
                })
                .expect("origins map is non-empty");

            let accepted = agg.n >= min_domains && agg.mean() >= min_score;
            OriginCall {
                read_id: read.id.clone(),
                origin: if accepted {
                    *origin
                } else {
                    Origin::Unclassified
                },
                score: agg.mean(),
                n_regions: agg.n,
                evalue: agg.best_evalue,
                env_from: if agg.env_from == usize::MAX {
                    0
                } else {
                    agg.env_from
                },
                env_to: agg.env_to,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slug_roundtrip() {
        for o in Origin::ALL {
            assert_eq!(Origin::from_slug(o.slug()), Some(o));
        }
    }

    fn reads(ids: &[&str]) -> Vec<IndexedRead> {
        ids.iter()
            .map(|id| IndexedRead {
                id: id.to_string(),
                sequence: "ACGTACGTACGT".to_string(),
            })
            .collect()
    }

    #[test]
    fn resolves_best_hit_per_read() {
        let hits = vec![
            HmmHit {
                read_id: "r1".into(),
                origin: Origin::Bacteria,
                region: "B01".into(),
                score: 100.0,
                evalue: 1e-20,
                env_from: 1,
                env_to: 1500,
            },
            HmmHit {
                read_id: "r1".into(),
                origin: Origin::Archaea,
                region: "A01".into(),
                score: 250.0,
                evalue: 1e-40,
                env_from: 1,
                env_to: 1500,
            },
        ];

        let calls = resolve(&reads(&["r1"]), &hits, 0.0, 1, RankMetric::Mean);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].origin, Origin::Archaea);
    }

    #[test]
    fn unclassifies_below_min_score() {
        let hits = vec![HmmHit {
            read_id: "r1".into(),
            origin: Origin::Bacteria,
            region: "B01".into(),
            score: 10.0,
            evalue: 0.5,
            env_from: 1,
            env_to: 100,
        }];

        let calls = resolve(&reads(&["r1"]), &hits, 50.0, 1, RankMetric::Mean);
        assert_eq!(calls[0].origin, Origin::Unclassified);
        // A weak hit keeps its score/coords for diagnostics.
        assert_eq!(calls[0].score, 10.0);
        assert_eq!(calls[0].env_to, 100);
    }

    #[test]
    fn no_hit_reads_become_unclassified() {
        // Two reads, only one has a hit; the other must still be reported.
        let hits = vec![HmmHit {
            read_id: "r1".into(),
            origin: Origin::Bacteria,
            region: "B01".into(),
            score: 100.0,
            evalue: 1e-20,
            env_from: 1,
            env_to: 1500,
        }];

        let calls = resolve(&reads(&["r1", "r2"]), &hits, 0.0, 1, RankMetric::Mean);
        assert_eq!(calls.len(), 2, "every input read must produce a call");
        // Order is preserved (input order).
        assert_eq!(calls[0].read_id, "r1");
        assert_eq!(calls[0].origin, Origin::Bacteria);
        assert_eq!(calls[0].n_regions, 1);
        assert_eq!(calls[1].read_id, "r2");
        assert_eq!(calls[1].origin, Origin::Unclassified);
        assert_eq!(calls[1].score, 0.0);
        assert_eq!(calls[1].n_regions, 0);
        assert_eq!(calls[1].env_from, 0);
        assert_eq!(calls[1].env_to, 0);
    }

    fn hit(origin: Origin, region: &str, score: f64) -> HmmHit {
        HmmHit {
            read_id: "r1".into(),
            origin,
            region: region.into(),
            score,
            evalue: 1e-10,
            env_from: 1,
            env_to: 1500,
        }
    }

    #[test]
    fn mean_score_beats_single_high_outlier() {
        // The real bacteria-vs-chloroplast failure mode: chloroplast has one
        // region that out-scores any single bacterial region, but bacteria wins
        // on mean across more regions. Old "max score wins" picked chloroplast;
        // mean-based ranking must pick bacteria.
        let hits = vec![
            hit(Origin::Bacteria, "B01", 46.0),
            hit(Origin::Bacteria, "B02", 45.0),
            hit(Origin::Bacteria, "B03", 44.0),
            hit(Origin::Chloroplast, "C01", 48.0), // higher single max
            hit(Origin::Chloroplast, "C02", 20.0),
        ];
        let calls = resolve(&reads(&["r1"]), &hits, 0.0, 2, RankMetric::Mean);
        assert_eq!(calls[0].origin, Origin::Bacteria);
        assert_eq!(calls[0].n_regions, 3);
        // envelope spans all matched bacterial regions
        assert_eq!(calls[0].env_from, 1);
        assert_eq!(calls[0].env_to, 1500);
    }

    #[test]
    fn min_domains_floor_suppresses_single_region_calls() {
        // One strong single-region hit, but min_domains = 2 → unclassified.
        let hits = vec![hit(Origin::Chloroplast, "C01", 90.0)];
        let calls = resolve(&reads(&["r1"]), &hits, 0.0, 2, RankMetric::Mean);
        assert_eq!(calls[0].origin, Origin::Unclassified);
        assert_eq!(calls[0].n_regions, 1);
        // and with the floor at 1 it is accepted
        let calls = resolve(&reads(&["r1"]), &hits, 0.0, 1, RankMetric::Mean);
        assert_eq!(calls[0].origin, Origin::Chloroplast);
    }

    #[test]
    fn rank_metric_changes_call_on_breadth_vs_mean_conflict() {
        // The real noisy-ONT shape: true bacteria has MORE regions at a slightly
        // lower mean; chloroplast has FEWER regions at a higher mean.
        //   B: 14 regions @ mean 34  (sum 476)
        //   C:  9 regions @ mean 36  (sum 324)
        let mut hits = Vec::new();
        for i in 0..14 {
            hits.push(hit(Origin::Bacteria, &format!("B{i:02}"), 34.0));
        }
        for i in 0..9 {
            hits.push(hit(Origin::Chloroplast, &format!("C{i:02}"), 36.0));
        }
        // Mean prefers the higher-mean organelle (the failure mode).
        let mean = resolve(&reads(&["r1"]), &hits, 0.0, 2, RankMetric::Mean);
        assert_eq!(mean[0].origin, Origin::Chloroplast);
        // Sum and Count both recover bacteria on breadth.
        let sum = resolve(&reads(&["r1"]), &hits, 0.0, 2, RankMetric::Sum);
        assert_eq!(sum[0].origin, Origin::Bacteria);
        let count = resolve(&reads(&["r1"]), &hits, 0.0, 2, RankMetric::Count);
        assert_eq!(count[0].origin, Origin::Bacteria);
    }
}
