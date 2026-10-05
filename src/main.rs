use anyhow::{Context, Result};
use clap::Parser;
use tracing::{info, Level};
use tracing_subscriber::EnvFilter;

use ssuplex::cli::Args;
use ssuplex::{extract, hmm, io, origin, summary, VERSION};

fn main() -> Result<()> {
    let args = Args::parse();
    init_logging(args.verbose);

    info!("SSUplex v{}", VERSION);
    info!("input: {}", args.input.display());
    info!("output prefix: {}", args.output.display());
    info!("threads: {}", args.threads);

    // Refuse to overwrite a previous run unless --force is given. The check
    // runs before any work; with --force, old files are removed only once the
    // new results are ready to be written.
    let previous: Vec<_> = summary::output_paths(&args)
        .into_iter()
        .filter(|p| p.exists())
        .collect();
    if !previous.is_empty() && !args.force {
        anyhow::bail!(
            "output files for prefix {} already exist (e.g. {}); use a new --output \
             prefix or pass --force to replace them",
            args.output.display(),
            previous[0].display()
        );
    }

    // Configure global rayon pool.
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("failed to initialise rayon thread pool")?;

    // Ensure the output directory exists before any writer runs.
    if let Some(dir) = args.output.parent() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir)
                .with_context(|| format!("failed to create output directory {}", dir.display()))?;
        }
    }

    // 1. Verify nhmmer is available and runnable.
    hmm::ensure_nhmmer_available()?;

    // 2. Discover HMM profiles for each origin.
    let profiles = hmm::load_origin_profiles(&args.hmm_dir).with_context(|| {
        format!(
            "failed to load HMM profiles from {}",
            args.hmm_dir.display()
        )
    })?;
    info!("loaded {} origin HMM profile(s)", profiles.len());

    // 3. Read all input sequences into memory, then scan them in chunks of
    //    --chunk-size against every origin profile.
    let reads = io::index_fasta(&args.input)
        .with_context(|| format!("failed to read {}", args.input.display()))?;
    info!("indexed {} input sequence(s)", reads.len());

    let hits = hmm::scan_all(&reads, &profiles, &args)?;

    if args.force {
        for p in summary::output_paths(&args).iter().filter(|p| p.exists()) {
            std::fs::remove_file(p)
                .with_context(|| format!("failed to remove previous output {}", p.display()))?;
        }
    }

    if args.debug_scores {
        summary::write_debug_scores(&args, &reads, &hits)?;
        summary::write_debug_regions(&args, &hits)?;
    }

    // 4. Resolve best origin per read.
    let calls = origin::resolve(&reads, &hits, args.min_score, args.min_domains, args.rank);
    let classified = calls
        .iter()
        .filter(|c| c.origin != origin::Origin::Unclassified)
        .count();
    info!(
        "classified {} of {} read(s) to an origin",
        classified,
        calls.len()
    );

    // 5. Extract SSU subsequences and write per-origin FASTAs + summary.
    let extracted = extract::extract_all(&reads, &calls)?;
    summary::write_all(&args, &reads, &calls, &extracted)?;

    info!("done.");
    Ok(())
}

fn init_logging(verbose: u8) {
    let level = match verbose {
        0 => Level::INFO,
        1 => Level::DEBUG,
        _ => Level::TRACE,
    };

    let filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(format!("ssuplex={level}")));

    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .init();
}
