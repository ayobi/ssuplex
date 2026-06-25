//! SSUplex: a fast Rust preprocessor for rRNA-marker eDNA workflows.
//!
//! Extracts SSU rRNA from environmental sequencing reads and sorts them by
//! origin (bacterial, archaeal, eukaryotic, mitochondrial, chloroplast) so
//! each can be routed to an appropriate downstream taxonomic classifier.
//! A faithful, parallelised Rust reimplementation of Metaxa2's SSU
//! extraction logic.
//!
//! The library is organised around a small pipeline:
//!
//! 1. [`io`]      — stream FASTA reads from disk in chunks
//! 2. [`hmm`]     — wrap `nhmmer` to scan each chunk against origin HMMs
//! 3. [`origin`]  — call best origin per read from collected hits
//! 4. [`extract`] — slice the SSU region out of each read using HMM coordinates
//! 5. [`summary`] — write origin-sorted FASTAs, coordinate TSV, and run summary
//!
//! The CLI in `main.rs` wires these together; library consumers can call them
//! individually for embedding into other pipelines.

pub mod cli;
pub mod extract;
pub mod hmm;
pub mod io;
pub mod origin;
pub mod summary;

/// Library version, surfaced in `--version` and the summary file.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
