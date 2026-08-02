extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::thread;

mod bamparser;
mod bedparser;
mod cli;

use crate::bamparser::BamParser;
use crate::bedparser::BedParser;
use crate::cli::ArgParser;

use crossbeam_channel::{select, unbounded, Receiver, Sender};
use indicatif::{ProgressBar, ProgressStyle};

type InputType = (String, u64, u64, Vec<(u64, u64)>);
type OutputType = String;

// Joins a slice of isize into a comma-separated string
fn join_isize(v: &[isize]) -> String {
    if v.is_empty() {
        return ".".to_string();
    }
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(",")
}

/// Chunk sorted BED-like regions by chromosome and proximity.
///
/// Entries are merged into the same chunk when they are on the same
/// chromosome and the gap between the current chunk's rightmost end
/// and the next entry's start is less than `chunksize`. Overlapping
/// entries (start <= current max_end) are always merged (distance 0).
///
/// Returns `(chrom, min_start, max_end, sub_intervals)` per chunk, where
/// `sub_intervals` holds the original `(start, end)` pairs in order.
fn chunk_regions(regions: Vec<(String, u64, u64)>, chunksize: u64) -> Vec<InputType> {
    let mut result: Vec<(String, u64, u64, Vec<(u64, u64)>)> = Vec::new();

    for (chrom, start, end) in regions {
        if let Some(last) = result.last_mut() {
            let same_chrom = last.0 == chrom;
            // saturating_sub so overlapping/contained intervals (start < max_end)
            // are treated as distance 0 instead of underflowing.
            let distance = start.saturating_sub(last.2);

            if same_chrom && distance < chunksize {
                last.2 = last.2.max(end);
                last.3.push((start, end));
                continue;
            }
        }

        // Start a new chunk.
        result.push((chrom, start, end, vec![(start, end)]));
    }

    result
}

fn main() -> std::io::Result<()> {
    let args = ArgParser::parse();
    let level = if args.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    pretty_env_logger::formatted_timed_builder()
        .filter_level(level)
        .init();

    args.validate();

    let mut file: Box<dyn Write> = match args.out {
        Some(ref path) => {
            let file = File::create(path).expect("Error creating output file");
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut m_parser = BedParser::new(&args.bed);

    let (sender, receiver): (Sender<Option<InputType>>, Receiver<Option<InputType>>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    info!("spawning {} threads", args.threads);
    for _ in 0..args.threads {
        let m_args = args.clone();
        let receiver = receiver.clone();
        let result_sender = result_sender.clone();
        thread::spawn(move || {
            let mut m_bam =
                BamParser::new(m_args.bam.clone(), m_args.reference.clone(), m_args.clone());
            for (chrom, start, end, sub_intervals) in receiver.into_iter().flatten() {
                let ret_str = if m_args.cov_only {
                    let data =
                        m_bam.extract_reads_coverage_only(&chrom, start, end, &sub_intervals);

                    data.iter()
                        .map(|(sub_start, sub_end, coverage)| {
                            format!("{}\t{}\t{}\t{}", chrom, sub_start, sub_end, coverage)
                        })
                        .collect::<Vec<_>>()
                        .join("\n")
                } else {
                    let data = m_bam.extract_reads_plup_fast(&chrom, start, end, &sub_intervals);

                    data.iter()
                        .map(|(sub_start, sub_end, coverage, deltas)| {
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                chrom,
                                sub_start,
                                sub_end,
                                coverage,
                                join_isize(&deltas[0]),
                                join_isize(&deltas[1]),
                                join_isize(&deltas[2]),
                            )
                        })
                        .collect::<Vec<_>>()
                        .join("\n")
                };
                result_sender.send(ret_str).unwrap();
            }
        });
    }

    let chunks = chunk_regions(m_parser.parse(), args.chunksize);
    let mut num_regions: u64 = 0;
    for region in chunks.into_iter() {
        sender.send(Some(region)).unwrap();
        num_regions += 1;
    }
    if num_regions == 0 {
        error!("No variants to be analyzed");
        std::process::exit(1);
    } else {
        info!("{} chunks to be analyzed", num_regions);
    }
    // Signal worker threads to exit
    for _ in 0..args.threads {
        sender.send(None).unwrap();
    }
    info!("collecting output");

    let sty =
        ProgressStyle::with_template(" [{elapsed_precise}] {bar:44.cyan/blue} > {pos} completed")
            .unwrap()
            .progress_chars("##-");
    let pbar = ProgressBar::new(num_regions);
    pbar.set_style(sty.clone());

    let mut collected: u64 = 0;
    // writeln!(file, "#chrom\tstart\tend\tcoverage\th0\th1\th2")?;
    loop {
        select! {
            recv(result_receiver) -> result => {
                match result {
                    Ok(ret_str) => {
                        writeln!(file, "{}", ret_str)?;
                        pbar.inc(1);
                        collected += 1;
                        if collected == num_regions {
                            break;
                        }
                    },
                    Err(e) => {
                        debug!("Problem {:?}", e);
                        break;
                    }
                }
            }
        }
    }
    pbar.finish();
    info!("finished");
    file.flush()
}
