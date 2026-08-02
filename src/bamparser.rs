use crate::cli::ArgParser;

use rust_htslib::bam::{self, record::Aux, IndexedReader, Read};

use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    args: ArgParser,
}

// Coverage, (hp0, hp1, hp2)
pub type DataTuple = (u64, u64, u64, [Vec<isize>; 3]);

impl BamParser {
    // Need the params to be passed here
    pub fn new(bam_name: PathBuf, reference: Option<PathBuf>, args: ArgParser) -> Self {
        let mut bam = IndexedReader::from_path(bam_name).unwrap();
        if let Some(ref_name) = reference {
            let _ = bam.set_reference(ref_name.clone());
        }
        BamParser { bam, args }
    }

    pub fn extract_reads_plup_fast(
            &mut self,
            chrom: &String,
            chunk_start: u64,
            chunk_end: u64,
            sub_intervals: &[(u64, u64)],
        ) -> Vec<DataTuple> {
            if let Err(e) = self.bam.fetch((&chrom, chunk_start, chunk_end)) {
                panic!(
                    "Unable to fetch bam {}:{}-{}\n{:?}",
                    chrom, chunk_start, chunk_end, e
                )
            };

        let anchor = self.args.anchor;
        let n = sub_intervals.len();
        let mut coverage: Vec<u64> = vec![0; n];
        let mut deltas: Vec<[Vec<isize>; 3]> =
            (0..n).map(|_| [Vec::new(), Vec::new(), Vec::new()]).collect();

        // Reused scratch buffers, reset per-read instead of reallocated.
        let mut diff_acc: Vec<isize> = vec![0; n];   // net indel delta within the core, deletions clipped to the core
        let mut left_gap: Vec<bool> = vec![false; n];  // a deletion/refskip touched the left flank => not continuous
        let mut right_gap: Vec<bool> = vec![false; n];
        let mut left_edit: Vec<u64> = vec![0; n];    // mismatch + insertion bases within the left flank
        let mut right_edit: Vec<u64> = vec![0; n];

        let mut alignment = bam::Record::new();

        // Nothing past this reference position can matter: the furthest any
        // sub-interval's right flank can reach is chunk_end + anchor.
        let scan_limit = chunk_end.saturating_add(self.args.anchor);

        while let Some(r) = self.bam.read(&mut alignment) {
            r.expect("Failed to parse record");

            if alignment.seq().is_empty()
                || alignment.is_unmapped()
                || alignment.mapq() < self.args.mapq
                || (alignment.flags() & self.args.mapflag) != 0
            {
                continue;
            }

            let ref_start = alignment.pos() as u64;

            // Only sub-intervals with >= `anchor` bp of overhang before their
            // left flank (ref_start <= sub_start - anchor) are candidates at all.
            // This is a cheap necessary-but-not-sufficient filter; the real
            // continuity/edit-distance check happens during the CIGAR walk below.
            let si_start = sub_intervals.partition_point(|&(s, _)| match s.checked_sub(self.args.anchor) {
                Some(threshold) => threshold < ref_start,
                None => true,
            });

            if si_start >= n {
                continue;
            }

            // Reset scratch state for this read.
            for j in si_start..n {
                diff_acc[j] = 0;
                left_gap[j] = false;
                right_gap[j] = false;
                left_edit[j] = 0;
                right_edit[j] = 0;
            }

            let mut read_pos = ref_start;
            let mut lo = si_start; // smallest index that might still need updates from later ops

            for &raw in alignment.raw_cigar() {
                let len = (raw >> 4) as u64;
                let opcode = raw & 0xf;

                // Op's reference-space span. Insertions are a zero-width
                // point event at the current reference position.
                let (op_start, op_end) = match opcode {
                    0 | 2 | 3 | 7 | 8 => (read_pos, read_pos + len), // M/D/N/=/X consume reference
                    _ => (read_pos, read_pos),                       // I/S/H/P etc: no ref span
                };

                // Permanently drop sub-intervals whose right flank ends at or
                // before this op starts -- nothing further in the read can
                // reach them, since ops only move rightward.
                while lo < n && sub_intervals[lo].1.saturating_add(anchor) <= op_start {
                    lo += 1;
                }

                for j in lo..n {
                    let (core_start, core_end) = sub_intervals[j];
                    let left_flank_start = core_start.saturating_sub(anchor);
                    let right_flank_end = core_end.saturating_add(anchor);

                    if left_flank_start >= op_end {
                        break; // this and all later sub-intervals start even further right
                    }

                    match opcode {
                        2 | 3 => {
                            // Deletion / RefSkip: a true gap in the read's sequence.
                            // Any overlap with a flank breaks "continuous" coverage
                            // of that flank outright.
                            if overlap_len(op_start, op_end, left_flank_start, core_start) > 0 {
                                left_gap[j] = true;
                            }
                            if overlap_len(op_start, op_end, core_end, right_flank_end) > 0 {
                                right_gap[j] = true;
                            }
                            // Only the portion of the deletion that actually falls
                            // inside the core counts toward the delta -- a deletion
                            // leading into/out of the region is clipped to its true
                            // overlap with the core.
                            let core_overlap = overlap_len(op_start, op_end, core_start, core_end);
                            if core_overlap > 0 {
                                diff_acc[j] -= core_overlap as isize;
                            }
                        }
                        1 => {
                            // Insertion: point event at `read_pos` (pre-advance).
                            if read_pos >= left_flank_start && read_pos < core_start {
                                left_edit[j] += len;
                            } else if read_pos >= core_end && read_pos < right_flank_end {
                                right_edit[j] += len;
                            } else if read_pos >= core_start && read_pos < core_end {
                                diff_acc[j] += len as isize;
                            }
                        }
                        8 => {
                            // Mismatch. Requires the aligner to emit extended
                            // CIGAR (=/X); plain 'M' can't be told apart from a
                            // true match without the reference sequence or MD tag.
                            left_edit[j] += overlap_len(op_start, op_end, left_flank_start, core_start);
                            right_edit[j] += overlap_len(op_start, op_end, core_end, right_flank_end);
                        }
                        _ => {}
                    }
                }

                read_pos = op_end;
                if read_pos > scan_limit {
                    break;
                }
            }

            let mut hp: Option<usize> = None;
            for j in si_start..n {
                let (_, core_end) = sub_intervals[j];
                let right_flank_end = core_end.saturating_add(anchor);

                if read_pos < right_flank_end {
                    continue; // didn't actually reach across the right flank
                }
                if left_gap[j] || right_gap[j] {
                    continue; // a deletion/refskip broke continuity over a flank
                }
                if self.args.anchor > 0 {
                    // edit_len / anchor < 0.10, done in integers as edit_len*10 < anchor
                    if left_edit[j] * self.args.max_edits >= self.args.anchor || right_edit[j] * self.args.max_edits >= self.args.anchor {
                        continue;
                    }
                }

                let hp_val = *hp.get_or_insert_with(|| Self::extract_hp(&alignment));
                coverage[j] += 1;
                deltas[j][hp_val].push(diff_acc[j]);
            }
        }

        sub_intervals
            .iter()
            .zip(coverage)
            .zip(deltas)
            .map(|(((s, e), cov), d)| (*s, *e, cov, d))
            .collect()
    }

    fn extract_hp(alignment: &bam::Record) -> usize {
        match alignment.aux(b"HP") {
            Ok(Aux::U8(value)) => value as usize,
            Ok(Aux::U16(value)) => value as usize,
            Ok(Aux::U32(value)) => value as usize,
            Ok(Aux::I32(value)) => value as usize,
            Ok(Aux::I16(value)) => value as usize,
            Ok(Aux::I8(value)) => value as usize,
            Ok(other) => panic!("Unexpected type for HP tag: {:?}", other),
            Err(_) => 0,
        }
    }
}

#[inline]
fn overlap_len(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> u64 {
    let lo = a_start.max(b_start);
    let hi = a_end.min(b_end);
    hi.saturating_sub(lo)
}
