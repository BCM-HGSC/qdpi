use crate::kmer::seq_to_kmer;
use crate::cli::ArgParser;
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{IndexedReader, Read},
    faidx,
};

use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    args: ArgParser,
}

impl BamParser {
    // Need the params to be passed here
    pub fn new(bam_name: PathBuf, ref_name: PathBuf,args: ArgParser) -> Self {
        let bam = IndexedReader::from_path(bam_name).unwrap();
        let reference = faidx::Reader::from_path(ref_name).unwrap();
        BamParser {
            bam,
            reference,
            args,
        }
    }
    
    // Return a vector of kfeats 
    pub fn extract_reads(&mut self, chrom: &String, start: u64, end: u64) -> Vec<Vec<f32>> {
        let ref_seq = remove_stretches(self.reference
                        .fetch_seq(chrom, start as usize, end as usize)
                        .unwrap().to_vec());

        let ref_kfeat = seq_to_kmer(&ref_seq, self.args.kmer, true);

        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        let start = start as usize;
        let end = end as usize;
        let mut alignments = Vec::new();

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

            //let read_name = String::from_utf8_lossy(alignment.qname()).to_string();
            // Do I need to utf8 this? I can probably keep it as bytes
            let read_sequence = alignment.seq().as_bytes();
            let cigar = alignment.cigar();

            let mut read_start: usize = 0;
            let mut read_len: usize = 0;
            let mut read_pos: usize = 0;
            let mut ref_pos: usize = alignment.pos() as usize;

            for op in &cigar {
                let m_len = op.len() as usize;
                match op.char() {
                    'M' | '=' | 'X' => {
                        if (ref_pos <= start) && (ref_pos + m_len) > start {
                            // We're over the start. Only happens once
                            let diff = start - ref_pos;
                            read_start = read_pos + diff;

                            // This match spans the whole region
                            if (ref_pos + m_len) >= end {
                                read_len += (ref_pos + m_len) - end;
                            } else {
                                // Or it stops inside
                                read_len += (ref_pos + m_len) - start;
                            }
                        } else if (ref_pos > start) && (ref_pos + m_len) < end {
                            // Match is within span
                            read_len += m_len;
                        } else if (ref_pos > start) && (ref_pos + m_len) >= end {
                            // Match is over the end
                            let trim = (ref_pos + m_len) - end;
                            read_len += m_len - trim;
                        }
                        ref_pos += m_len;
                        read_pos += m_len;
                    }
                    'I' => {
                        if (ref_pos >= start) && (ref_pos <= end) {
                            read_len += m_len;
                        }
                        read_pos += m_len;
                    },
                    'D' => {
                        // If this deletion spans the start, we need to set the start
                        if (ref_pos <= start) && (ref_pos + m_len) >= start {
                            let diff = start - ref_pos;
                            read_start = read_pos + diff + 1; // the next base will be the start
                        }
                        ref_pos += m_len;
                    }
                    'S' => {
                        // Skip the soft-clipped bases
                        read_pos += m_len;
                    }
                    _ => {}
                }
                // Don't need to keep trimming
                if ref_pos > end {
                    break;
                }
            }
            
            let seq = if read_len != 0 {
                let m_seq = remove_stretches(read_sequence[read_start..(read_start + read_len)].to_vec());
                m_seq
            } else {
                vec![]
            };
            debug!("{:?}", String::from_utf8_lossy(&seq).to_string());
            let mut m_kfeat = seq_to_kmer(&seq, self.args.kmer, false);
            m_kfeat
            .iter_mut()
            .zip(ref_kfeat.iter())
            .for_each(|(x, y)| *x += y);


            alignments.push(m_kfeat);
        }

        alignments
    }
}


fn remove_stretches(vector: Vec<u8>) -> Vec<u8> {
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

