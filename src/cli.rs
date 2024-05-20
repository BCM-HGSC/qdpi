extern crate pretty_env_logger;

use clap::Parser;

#[derive(Parser, Clone, Debug)]
#[command(author = "ACEnglish", version)]
pub struct ArgParser {
    /// Reads to analyze
    #[arg(short, long)]
    pub bam: std::path::PathBuf,

    /// Regions to analyze
    #[arg(long)]
    pub bed: std::path::PathBuf,

    /// Output tsv
    #[arg(short, long)]
    pub out: std::path::PathBuf,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Alignments with flag matching this value are ignored
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,

    /// Number of threads
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    /// Verbose logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,

    /// Buffer around regions to plup
    #[arg(long, default_value_t = 50)]
    pub buffer: u64,
}

impl ArgParser {
    /// Validate command line arguments
    pub fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.bam.exists() {
            error!("--bam does not exist");
            is_ok = false;
        } else if !self.bam.is_file() {
            error!("--bam is not a file");
            is_ok = false;
        }

        if !self.bed.exists() {
            error!("--bed does not exist");
            is_ok = false;
        } else if !self.bed.is_file() {
            error!("--bed is not a file");
            is_ok = false;
        }

        is_ok
    }
}
