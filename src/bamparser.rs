use crate::cli::ArgParser;
use rust_htslib::bam::{self, record::Cigar};
use rust_htslib::bam::record::CigarStringView;

use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::pileup::Indel,
    bam::{IndexedReader, Read},
};

use std::collections::HashMap;
use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    args: ArgParser,
}

type DataTuple = (u64, Vec<(u64, Vec<String>)>, Vec<isize>);
type DataTuple2 = (u64, Vec<(u64, Vec<isize>)>, Vec<isize>);
impl BamParser {
    // Need the params to be passed here
    pub fn new(bam_name: PathBuf, args: ArgParser) -> Self {
        let bam = IndexedReader::from_path(bam_name).unwrap();
        BamParser { bam, args }
    }

    pub fn extract_reads_plup(&mut self, chrom: &String, start: u64, end: u64) -> DataTuple {
        let mut tot_cov: u64 = 0;
        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        // track the changes made by each read
        let mut m_reads = HashMap::<Vec<u8>, isize>::new();
        let mut plups = HashMap::<u64, Vec<String>>::new();
        for pileup in self.bam.pileup().flatten() {
            let m_pos: u64 = pileup.pos().into();
            // We got to truncate the pileup
            if m_pos < start {
                continue;
            }
            // Do this separately so we don't waste time downstream
            if end < m_pos {
                break;
            }
            for alignment in pileup.alignments() {
                if alignment.record().seq().is_empty()
                    || alignment.qpos().is_none()
                    || alignment.record().mapq() < self.args.mapq
                    || (alignment.record().flags() & self.args.mapflag) != 0
                    || (!((alignment.record().reference_start() as u64) < start
                        && (alignment.record().reference_end() as u64) > end))
                {
                    continue;
                }

                tot_cov += 1;
                let (m_var, m_seq) = match alignment.indel() {
                    Indel::Ins(size) => {
                        let qpos = alignment.qpos().unwrap();
                        let seq = alignment.record().seq().as_bytes()[qpos..(qpos + size as usize)]
                            .to_vec();
                        let seq = std::str::from_utf8(&seq).unwrap_or("").to_string();
                        (size as isize, seq)
                    }
                    Indel::Del(size) => (-(size as isize), format!("-{}", size)),
                    _ => continue,
                };
                plups.entry(m_pos - start).or_default().push(m_seq);
                let qname = alignment.record().qname().to_owned();
                *m_reads.entry(qname).or_insert(0) += m_var;
            }
        }

        let coverage = (tot_cov / (end - start)).max(m_reads.len() as u64);
        let plups = plups.into_iter().collect();
        let deltas: Vec<_> = m_reads.into_values().collect();

        (coverage, plups, deltas)
    }

    pub fn extract_reads_plup_fast(&mut self, chrom: &String, start: u64, end: u64) -> DataTuple2 {
        let mut tot_cov: u64 = 0;
        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        // track the changes made by each read
        let mut deltas = Vec::new();
        let mut plups = HashMap::<u64, Vec<isize>>::new();
        let mut coverage = 0;
        for alignment in self.bam.records() {
            let alignment = alignment.expect("Right?");
            if alignment.seq().is_empty()
                || alignment.is_unmapped()
                || alignment.mapq() < self.args.mapq
                || (alignment.flags() & self.args.mapflag) != 0
                || (!((alignment.reference_start() as u64) < start
                    && (alignment.reference_end() as u64) > end))
            {
                continue;
            }
            coverage += 1;
            let read_start = alignment.pos() as u64;
            let cigar = alignment.cigar();
            let mut read_pos = read_start;
            let mut total_diff: isize = 0;

            for cig in cigar.iter() {
                match cig {
                    Cigar::Ins(len) => {
                        if read_pos >= start && read_pos <= end {
                            let len = len.clone() as isize;
                            total_diff += len;
                            plups.entry(read_pos - start).or_default().push(len);
                        }
                    }
                    Cigar::Del(len) => {
                        let len = len.clone() as isize;
                        if read_pos >= start && read_pos <= end {
                            total_diff -= len;
                            plups.entry(read_pos - start).or_default().push(-len);
                        }
                        read_pos += len as u64;
                    }
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        read_pos += *len as u64;
                    }
                    Cigar::SoftClip(len)
                    | Cigar::HardClip(len)
                    | Cigar::Pad(len)
                    | Cigar::RefSkip(len) => {
                        // These do not affect read position for the query sequence
                    }
                    _ => {}
                }
            }


            deltas.push(total_diff);
        }

        let plups = plups.into_iter().collect();

        (coverage, plups, deltas)
    }
}
