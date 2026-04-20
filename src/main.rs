use anyhow::Result;
use clap::Parser;
use log::info;

use markdup_wea::barcode_tags::BarcodeTags;
use markdup_wea::markdup;

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

    /// BARCODE_TAG — 2-char SAM aux tag carrying the whole-pair UMI/barcode
    /// (Picard `BARCODE_TAG`). Typical value: RX. When set, reads at the same
    /// coordinates with different barcode values are NOT grouped as duplicates.
    #[arg(long, value_parser = parse_aux_tag)]
    barcode_tag: Option<String>,

    /// READ_ONE_BARCODE_TAG — barcode tag on the firstOfPair record
    /// (Picard `READ_ONE_BARCODE_TAG`). 2-char SAM aux tag.
    #[arg(long, value_parser = parse_aux_tag)]
    read_one_barcode_tag: Option<String>,

    /// READ_TWO_BARCODE_TAG — barcode tag on the !firstOfPair record
    /// (Picard `READ_TWO_BARCODE_TAG`). 2-char SAM aux tag.
    #[arg(long, value_parser = parse_aux_tag)]
    read_two_barcode_tag: Option<String>,
}

/// Validate a SAM aux tag identifier. Per SAMv1 §1.5, tags match
/// `[A-Za-z][A-Za-z0-9]`. Accept only that pattern (e.g. RX, BC, B2).
fn parse_aux_tag(s: &str) -> Result<String, String> {
    let b = s.as_bytes();
    if b.len() != 2 {
        return Err(format!("SAM aux tag must be 2 chars, got {:?}", s));
    }
    if !b[0].is_ascii_alphabetic() {
        return Err(format!(
            "SAM aux tag first char must be a letter, got {:?}",
            s
        ));
    }
    if !b[1].is_ascii_alphanumeric() {
        return Err(format!(
            "SAM aux tag second char must be a letter or digit, got {:?}",
            s
        ));
    }
    Ok(s.to_string())
}

fn to_tag_bytes(s: &str) -> [u8; 2] {
    let b = s.as_bytes();
    [b[0], b[1]]
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

    // Own the 2-byte tag arrays for the whole run; `BarcodeTags` borrows them.
    let barcode_bytes = cli.barcode_tag.as_deref().map(to_tag_bytes);
    let read_one_bytes = cli.read_one_barcode_tag.as_deref().map(to_tag_bytes);
    let read_two_bytes = cli.read_two_barcode_tag.as_deref().map(to_tag_bytes);
    let barcode_tags = BarcodeTags {
        barcode: barcode_bytes.as_ref(),
        read_one: read_one_bytes.as_ref(),
        read_two: read_two_bytes.as_ref(),
    };

    markdup::run(
        &cli.input,
        cli.output.as_deref(),
        cli.metrics.as_deref(),
        cli.threads,
        cli.remove_duplicates,
        cli.assume_sort_order.as_deref(),
        barcode_tags,
    )
}
