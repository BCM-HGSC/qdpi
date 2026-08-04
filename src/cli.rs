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

    /// Reference path (for crams)
    #[arg(long)]
    pub reference: Option<std::path::PathBuf>,

    /// Output bed-like (default: stdout)
    #[arg(short, long)]
    pub out: Option<std::path::PathBuf>,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 60)]
    pub mapq: u8,

    /// Alignments with flag matching this value are ignored
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,

    /// Alignments must span up/down stream Nbp flanks around region
    #[arg(long, default_value_t = 200)]
    pub flank: u64,

    /// Filter alignments that don't match to N% of flanks
    #[arg(long, default_value_t = 0.90)]
    pub min_anchor: f32,

    /// Number of threads
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    /// Verbose logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,

    /// Chunksize for optimized IO - set to a little longer than read-length
    #[arg(long, default_value_t = 20000)]
    pub chunksize: u64,
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

        if let Some(refname) = &self.reference {
            if !refname.exists() {
                error!("--reference does not exist");
                is_ok = false;
            } else if !refname.is_file() {
                error!("--reference is not a file");
                is_ok = false;
            }
        }

        is_ok
    }
}
