use anyhow::Result;
use clap::Parser;
use log::info;

mod dupset;
mod groups;
mod io;
mod markdup;
mod metrics;
mod pending_mates;
mod position;
mod scan;
mod scoring;

#[derive(Parser)]
#[command(name = "markdup-wea")]
#[command(about = "Fast BAM duplicate marker — Picard-compatible, Rust-powered")]
#[command(version)]
struct Cli {
    /// Input BAM file (coordinate-sorted). Use - for stdin.
    input: String,

    /// Output BAM file [default: stdout]
    #[arg(short, long)]
    output: Option<String>,

    /// Picard-format metrics file
    #[arg(short = 'M', long)]
    metrics: Option<String>,

    /// I/O threads for BAM reading/writing
    #[arg(short = '@', long = "threads", default_value = "1")]
    threads: u32,

    /// Exclude duplicates from output instead of flagging
    #[arg(long)]
    remove_duplicates: bool,

    /// Override header sort order check
    #[arg(long)]
    assume_sort_order: Option<String>,
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    info!(
        "markdup-wea v{} — input: {}, threads: {}",
        env!("CARGO_PKG_VERSION"),
        cli.input,
        cli.threads
    );

    markdup::run(
        &cli.input,
        cli.output.as_deref(),
        cli.metrics.as_deref(),
        cli.threads,
        cli.remove_duplicates,
        cli.assume_sort_order.as_deref(),
    )
}
