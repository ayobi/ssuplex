//! Integration tests.
//!
//! These run without HMMER. End-to-end runs (detection, clipping, overwrite
//! handling with real outputs) are covered by `dev/smoke.sh`.

use assert_cmd::Command;
use predicates::prelude::*;

#[test]
fn shows_help() {
    let mut cmd = Command::cargo_bin("ssuplex").unwrap();
    cmd.arg("--help").assert().success();
}

#[test]
fn shows_version() {
    let mut cmd = Command::cargo_bin("ssuplex").unwrap();
    cmd.arg("--version").assert().success();
}

#[test]
fn refuses_to_overwrite_previous_outputs() {
    let dir = tempfile::tempdir().unwrap();
    let previous = dir.path().join("sample.summary.txt");
    std::fs::write(&previous, "previous run").unwrap();
    Command::cargo_bin("ssuplex")
        .unwrap()
        .args(["-i", "missing.fasta", "-o"])
        .arg(dir.path().join("sample"))
        .assert()
        .failure()
        .stderr(predicate::str::contains("--force"));
    assert_eq!(std::fs::read_to_string(&previous).unwrap(), "previous run");
}

#[test]
fn force_passes_the_guard_without_deleting_before_a_successful_scan() {
    // With --force the guard is skipped; the run then fails on the missing
    // input, and the previous results must still be there.
    let dir = tempfile::tempdir().unwrap();
    let previous = dir.path().join("sample.summary.txt");
    std::fs::write(&previous, "previous run").unwrap();
    Command::cargo_bin("ssuplex")
        .unwrap()
        .args(["-i", "missing.fasta", "--force", "-o"])
        .arg(dir.path().join("sample"))
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exist").not());
    assert!(previous.exists());
}
