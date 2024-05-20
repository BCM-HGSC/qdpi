use crate::cli::ArgParser;
use crate::kmer::seq_to_kmer;
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


    // Return a vector of kfeats
    pub fn extract_reads(&mut self, chrom: &String, start: u64, end: u64) -> Locus {
        let ref_seq = remove_homopolymer(
            self.reference
                .fetch_seq(chrom, start as usize, end as usize)
                .unwrap()
                .to_vec(),
        );

        let ref_kfeat = seq_to_kmer(&ref_seq, self.args.kmer, true);

        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        let mut result = Locus::new(chrom, start, end);
        let start = start as usize;
        let end = end as usize;

        // I think I can optimize this to read into a single object or something?
        for alignment in self.bam.records() {
            let alignment = alignment.expect("Failed to read alignment");
            // Skip records without sequence, below min mapq, matching the flag,
            // or partially into our window. Will revisit partial when I can turn a Haplotype into a
            // single-path graph
            if alignment.seq().is_empty()
                || alignment.mapq() < self.args.mapq
                || (alignment.flags() & self.args.mapflag) != 0
                || !((alignment.reference_start() as usize) < start
                    && (alignment.reference_end() as usize) > end)
            {
                continue;
            }

            // Do I need to utf8 this? I can probably keep it as bytes
            let read_sequence = alignment.seq().as_bytes();
            let cigar = alignment.cigar();

            let mut read_start: usize = 0;
            let mut read_len: usize = 0;
            let mut read_pos: usize = 0;
            let mut ref_pos: usize = alignment.pos() as usize;

            for op in &cigar {
                if ref_pos >= end {
                    break;
                }

                let op_len = op.len() as usize;
                match op.char() {
                    'M' | '=' | 'X' => {
                        if (ref_pos <= start) && (ref_pos + op_len) > start {
                            // We're over the start. Only happens once
                            let diff = start - ref_pos;
                            read_start = read_pos + diff;

                            // This match spans the whole region
                            if (ref_pos + op_len) >= end {
                                read_len += end - start;
                            } else {
                                // Or it stops inside
                                read_len += (ref_pos + op_len) - start;
                            }
                        } else if (ref_pos > start) && (ref_pos + op_len) < end {
                            // Match is within span
                            read_len += op_len;
                        } else if (ref_pos > start) && (ref_pos + op_len) >= end {
                            // Match is over the end
                            let trim = (ref_pos + op_len) - end;
                            read_len += op_len - trim;
                        }
                        ref_pos += op_len;
                        read_pos += op_len;
                    }
                    'D' | 'N' => {
                        // The rest is deleted
                        if (ref_pos + op_len) > end {
                            break;
                        }
                        if (ref_pos <= start) && (ref_pos + op_len) >= start {
                            let diff = start - ref_pos;
                            read_start = read_pos + diff + 1; // the next base will be the start
                        }
                        ref_pos += op_len;
                    }
                    'I' | 'S' | 'H' => {
                        if (ref_pos >= start) && (ref_pos <= end) {
                            read_len += op_len;
                        }
                        read_pos += op_len;
                    }
                    _ => {}
                }
            }

            // I have some kind of trimming error because I get start+len > read sometimes
            if (read_start + read_len) >= read_sequence.len() {
                continue;
            }

            result.n_reads += 1;
            if (read_len as isize).abs_diff((end - start) as isize) < 20 {
                continue; // must be at least 20bp diff
            }

            if read_len == 0 {
                continue
            }

            let seq = remove_homopolymer(read_sequence[read_start..(read_start + read_len)].to_vec());
            //let read_name = String::from_utf8_lossy(alignment.qname()).to_string();
            //debug!("{:?} {} {} {} {:?}", read_name, read_sequence.len(), alignment.pos(),
                    //seq.len(), String::from_utf8_lossy(&seq).to_string());
            let mut m_kfeat = seq_to_kmer(&seq, self.args.kmer, false);
            m_kfeat
                .iter_mut()
                .zip(ref_kfeat.iter())
                .for_each(|(x, y)| *x += y);

            //*result.read_lengths.entry(read_len).or_insert(0) += 1;
            // Debug 
            result.read_seqs.push(seq);
            result.read_lengths.push(read_len);
            result.kfeats.push(m_kfeat);
        }

        result
    }
}

fn remove_homopolymer(vector: Vec<u8>) -> Vec<u8> {
    let maxspan = 5;
    let mut result = Vec::new();
    let mut count = 0;
    let mut prev_byte = None;

    for byte in vector {
        match prev_byte {
            Some(prev) if prev == byte => {
                count += 1;
                if count < maxspan {
                    result.push(byte);
                }
            }
            _ => {
                count = 1;
                result.push(byte);
            }
        }
        prev_byte = Some(byte);
    }

    result
}
