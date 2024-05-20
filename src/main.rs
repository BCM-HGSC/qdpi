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

type InputType = (String, u64, u64);
type OutputType = (String, u64, u64, String);

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

    let fh = File::create(args.out.clone())?;
    let mut file = BufWriter::new(fh);

    let mut m_parser = BedParser::new(&args.bed);

    let (sender, receiver): (Sender<Option<InputType>>, Receiver<Option<InputType>>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    info!("spawning {} threads", args.threads);
    for _ in 0..args.threads {
        let m_args = args.clone();
        let receiver = receiver.clone();
        let result_sender = result_sender.clone();
        thread::spawn(move || {
            let mut m_bam = BamParser::new(m_args.bam.clone(), m_args.clone());
            for (chrom, mut start, mut end) in receiver.into_iter().flatten() {
                start -= m_args.buffer;
                end += m_args.buffer;
                let (reads, plups, coverage) = m_bam.extract_reads_plup(&chrom, start, end);
                let deltas: Vec<_> = reads.values().collect();
                let mut alt_ps = Vec::<(usize, Vec<String>)>::new();
                for (p, _, a) in plups {
                    let shift_p = (p - start) as usize;
                    if !a.is_empty() {
                        alt_ps.push((shift_p, a))
                    }
                }
                let data = (coverage, alt_ps, deltas);
                let json_str = match serde_json::to_string(&data) {
                    Ok(json) => json,
                    Err(e) => {
                        eprintln!("Error serializing to JSON: {}", e);
                        "".to_string()
                    }
                };

                result_sender.send((chrom, start, end, json_str)).unwrap();
            }
        });
    }

    let mut num_regions: u64 = 0;
    for region in m_parser.parse().into_iter() {
        sender.send(Some(region)).unwrap();
        num_regions += 1;
    }
    if num_regions == 0 {
        error!("No variants to be analyzed");
        std::process::exit(1);
    } else {
        info!("{} regions to be analyzed", num_regions);
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
    loop {
        select! {
            recv(result_receiver) -> result => {
                match result {
                    Ok((chrom, start, end, json_str)) => {
                        write!(file, "{}\t{}\t{}\t{}\n", chrom, start, end, json_str)?;
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
