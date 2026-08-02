use crate::cli::ArgParser;

use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{IndexedReader, Read},
    bam::record::{Aux, Cigar}
};

use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    args: ArgParser,
}

// Coverage, (hp0, hp1, hp2)
type DataTuple = (u64, [Vec<isize>; 3]);

impl BamParser {
    // Need the params to be passed here
    pub fn new(bam_name: PathBuf, reference: Option<PathBuf>, args: ArgParser) -> Self {
        let mut bam = IndexedReader::from_path(bam_name).unwrap();
        if let Some(ref_name) = reference {
            let _ = bam.set_reference(ref_name.clone());
        }
        BamParser { bam, args }
    }
    
    pub fn extract_reads_plup_fast(&mut self, chrom: &String, start: u64, end: u64) -> DataTuple {
        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        let mut coverage = 0;
        let mut deltas = [Vec::new(), Vec::new(), Vec::new()];

        for alignment in self.bam.records() {
            let alignment = alignment.expect("Right?");

            if alignment.seq().is_empty()
                || alignment.is_unmapped()
                || alignment.mapq() < self.args.mapq
                || (alignment.flags() & self.args.mapflag) != 0
            {
                continue;
            }

            let ref_start = alignment.pos() as u64;
            if ref_start >= start {
                continue; // can't possibly span the region from the left
            }

            // Single pass: determines span AND accumulates the indel delta.
            let mut read_pos = ref_start;
            let mut total_diff: isize = 0;
            let mut spans_region = false;

            for &raw in alignment.raw_cigar() {
                let len = (raw >> 4) as u64;
                match raw & 0xf {
                    1 => { // Ins
                        if read_pos >= start && read_pos <= end {
                            total_diff += len as isize;
                        }
                    }
                    2 => { // Del
                        if read_pos >= start && read_pos <= end {
                            total_diff -= len as isize;
                        }
                        read_pos += len;
                    }
                    0 | 7 | 8 => { // Match/Equal/Diff
                        read_pos += len;
                    }
                    _ => {}
                }
                if read_pos > end {
                    spans_region = true;
                    break;
                }
            }

            if !spans_region {
                continue;
            }
            coverage += 1;

            let hp = match alignment.aux(b"HP") {
                Ok(Aux::U8(value)) => value as usize,
                Ok(Aux::U16(value)) => value as usize,
                Ok(Aux::U32(value)) => value as usize,
                Ok(Aux::I32(value)) => value as usize,
                Ok(Aux::I16(value)) => value as usize,
                Ok(Aux::I8(value)) => value as usize,
                Ok(other) => panic!("Unexpected type for HP tag: {:?}", other),
                Err(_) => 0,
            };

            deltas[hp].push(total_diff);
        }

        (coverage, deltas)
    }

    pub fn __extract_reads_plup_fast(&mut self, chrom: &String, start: u64, end: u64) -> DataTuple {
        if let Err(e) = self.bam.fetch((&chrom, start, end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        let mut coverage = 0;
        // track the changes made by each read
        // By HP (0, 1, 2)
        let mut deltas = [Vec::new(), Vec::new(), Vec::new()];

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
            let cigar = alignment.cigar();
            let mut read_pos = alignment.pos() as u64;
            let mut total_diff: isize = 0;

            for cig in cigar.iter() {
                match cig {
                    Cigar::Ins(len) => {
                        if read_pos >= start && read_pos <= end {
                            total_diff += *len as isize;
                        }
                    }
                    Cigar::Del(len) => {
                        if read_pos >= start && read_pos <= end {
                            total_diff -= *len as isize;
                        }
                        read_pos += *len as u64;
                    }
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        read_pos += *len as u64;
                    }
                    _ => {}
                }

                if read_pos > end {
                    break;
                }
            }

            let hp = match alignment.aux(b"HP") {
                Ok(Aux::U8(value)) => value as usize,
                Ok(Aux::U16(value)) => value as usize,
                Ok(Aux::U32(value)) => value as usize,
                Ok(Aux::I32(value)) => value as usize,
                Ok(Aux::I16(value)) => value as usize,
                Ok(Aux::I8(value)) => value as usize,
                Ok(other) => {
                    panic!("Unexpected type for HP tag: {:?}", other);
                }
                Err(_) => 0,
            };


            deltas[hp].push(total_diff);
        }

        (coverage, deltas)
    }
}
