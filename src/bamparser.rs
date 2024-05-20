use crate::cli::ArgParser;
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::pileup::Indel,
    bam::{IndexedReader, Read},
    faidx,
};

use crate::locus::Locus;
use std::path::PathBuf;
use std::collections::HashMap;

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    args: ArgParser,
}

impl BamParser {
    // Need the params to be passed here
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, args: ArgParser) -> Self {
        let bam = IndexedReader::from_path(bam_name).unwrap();
        let reference = faidx::Reader::from_path(ref_name).unwrap();
        BamParser {
            bam,
            reference,
            args,
        }
    }

    pub fn extract_reads_plup(&mut self, chrom: &String, start: u64, end: u64) -> (HashMap::<String, isize>, Vec<(u64, u32, Vec<String>)>, u64) {
        let mut tot_cov: u64 = 0;
        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        // track the changes made by each read
        let mut m_reads = HashMap::<String, isize>::new();
        let mut plups = Vec::<(u64, u32, Vec<String>)>::new();
        for pileup in self.bam.pileup() {
            if pileup.is_err() {
                continue;
            }
            let pileup = pileup.unwrap();
            let m_pos: u64 = pileup.pos().into();
            // We got to truncate the pileup
            if m_pos < start {
                continue;
            }
            // Do this separately so we don't waste time downstream
            if end < m_pos {
                break;
            }
            let mut ps = Vec::new();
            for alignment in pileup.alignments() {
                if alignment.record().seq().is_empty() || alignment.qpos().is_none()
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
                        let seq =
                            alignment.record().seq().as_bytes()[qpos..(qpos + size as usize)].to_vec();
                        let seq = std::str::from_utf8(&seq).unwrap_or("").to_string();
                        (size as isize, seq)
                    },
                    Indel::Del(size) => {
                        (-(size as isize), format!("-{}", size))

                    }
                    _ => continue,
                };
                ps.push(m_seq);
                let qname = String::from_utf8_lossy(alignment.record().qname()).to_string();
                //let qname = alignment.record().qname().to_owned();
                *m_reads.entry(qname).or_insert(0) += m_var;
            }
            plups.push((m_pos, pileup.depth(), ps));
        }

        let coverage = (tot_cov / (end - start)).max(m_reads.len() as u64);
        (m_reads, plups, coverage)
    }
}
