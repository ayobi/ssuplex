//! Integration tests.
//!
//! Most of these will require HMMER on `$PATH` and small bundled HMM/FASTA
//! fixtures. Add fixtures under `tests/data/` as the test suite grows.

use assert_cmd::Command;

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
