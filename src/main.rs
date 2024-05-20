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
mod kmeans;
mod kmer;
mod locus;

use crate::bamparser::BamParser;
use crate::bedparser::BedParser;
use crate::cli::ArgParser;
use crate::kmeans::kmeans;
use crate::kmeans::Cluster;

use serde_json::json;
use crossbeam_channel::{select, unbounded, Receiver, Sender};
use indicatif::{ProgressBar, ProgressStyle};

type InputType = (String, u64, u64);
type OutputType = (String, u64, u64, String);

fn calinski_harabasz_index(clusters: &[Cluster]) -> f32 {
    let n_clusters = clusters.len();
    let n_points = clusters
        .iter()
        .map(|cluster| cluster.points.len())
        .sum::<usize>();
    let n_features = clusters[0].centroid.len();
    if n_points <= n_clusters {
        return f32::NAN;
    }
    // Calculate centroids of centroids
    let centroid_of_centroids: Vec<f32> = clusters
        .iter()
        .flat_map(|cluster| cluster.centroid.iter())
        .fold(vec![0.0; n_features], |mut acc, &x| {
            acc.iter_mut().for_each(|a| *a += x);
            acc
        })
        .iter()
        .map(|&sum| sum / n_features as f32)
        .collect();

    // Calculate between-cluster dispersion
    let between_cluster_dispersion: f32 =
        clusters
            .iter()
            .filter(|v| !v.points.is_empty())
            .fold(0.0, |acc, cluster| {
                let cluster_size = cluster.points.len() as f32;
                acc + cluster_size * squared_distance(&cluster.centroid, &centroid_of_centroids)
            });

    // Calculate within-cluster dispersion
    let within_cluster_dispersion: f32 =
        clusters
            .iter()
            .filter(|v| !v.points.is_empty())
            .fold(0.0, |acc, cluster| {
                acc + cluster.points.iter().fold(0.0, |acc, point| {
                    acc + squared_distance(point, &cluster.centroid)
                })
            });

    // Within is perfect, so just return between dispersion?
    if within_cluster_dispersion == 0.0 {
        between_cluster_dispersion / (n_clusters - 1) as f32
    } else {
        // Calculate Calinski-Harabasz index
        (between_cluster_dispersion / (n_clusters - 1) as f32)
            / (within_cluster_dispersion / (n_points - n_clusters) as f32)
    }
}

fn squared_distance(point1: &[f32], point2: &[f32]) -> f32 {
    point1
        .iter()
        .zip(point2.iter())
        .map(|(&x, &y)| (x - y).powi(2))
        .sum()
}

/// Canberra distance of featurized kmers
pub fn seqsim(a: &[f32], b: &[f32], mink: f32) -> f32 {
    let mut deno: f32 = 0.0;
    let mut neum: f32 = 0.0;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let d = x.abs() + y.abs();
        if d <= mink {
            continue;
        }
        deno += d;

        neum += (x - y).abs();
    }

    // no kmers
    if deno == 0.0 {
        return 0.0;
    }

    // identical
    if neum == 0.0 {
        return 1.0;
    }

    1.0 - (neum / deno)
}

fn __sum_of_squared_distances(centroid: &[f32], values: &[Vec<f32>]) -> f32 {
    let mut sum = 0.0;

    for value in values {
        sum += seqsim(centroid, value, 1.0);
    }

    sum / values.len() as f32
}

fn __find_elbow_point(distances: &[f32]) -> usize {
    let mut best_val = distances[0];
    let mut elbow_point = 1;

    for i in distances.iter().skip(1) {
        if *i <= best_val {
            break;
        }
        if *i > best_val {
            best_val = *i;
            elbow_point += 1;
        }
    }

    elbow_point
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
            let mut m_bam = BamParser::new(m_args.bam.clone(), m_args.reference.clone(), m_args.clone());
            for (chrom, mut start, mut end) in receiver.into_iter().flatten() {
                start -= m_args.buffer;
                end += m_args.buffer;
                let (reads, plups, coverage) = m_bam.extract_reads_plup(&chrom, start, end);
                let deltas:Vec<_> = reads.values().collect();
                //let mut pos_cov = vec![0; (end - start + 1) as usize];
                let mut alt_ps = Vec::<(usize, Vec<String>)>::new();
                for (p, c, a) in plups {
                    let shift_p = (p - start) as usize;
                    //pos_cov[shift_p] += c;
                    if !a.is_empty() {
                        alt_ps.push((shift_p, a))
                    }
                }
                let data = (coverage, alt_ps, deltas);
                let json_str = match serde_json::to_string(&data) {
                    Ok(json) => json,
                    Err(e) => { eprintln!("Error serializing to JSON: {}", e); "".to_string() },
                };
                /*let mut json_obj: serde_json::Map<String, serde_json::Value> = result
                    .iter()
                    .map(|(k, v)| (String::from_utf8_lossy(k).to_string(), json!(v)))
                    .collect();
                json_obj["P"] = plups;*/

                /*let json_str = serde_json::to_string(&json_obj).unwrap();*/

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
    }
     else {
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

/*
 *    for (chrom, start, end) in m_parser.parse().into_iter() {
        info!("{}:{}-{}", chrom, start, end);

        let mut result = m_bam.extract_reads(&chrom, start, end);
        if result.kfeats.is_empty() {
            continue;
        }

        let mut best_score = 0.0;
        let mut best_k = 1;
        let mut points = vec![];
        for k in 2..10 {
            let clusts = kmeans(&result.kfeats, k);
            let idx = calinski_harabasz_index(&clusts);
            debug!("{} {}", k, idx);
            if !idx.is_finite() || idx < best_score {
                break;
            }
            best_score = idx;
            best_k = k;
            points = clusts;
        }
        result.n_clusters = best_k;
        result.cluster_score = best_score;
        result.clusters = points.iter().map(|p| p.points_idx.clone()).collect();

        result.make_data(&mut file)?;
        //result.make_verbose_data(&mut file)?;
    }

    file.flush()
*/
