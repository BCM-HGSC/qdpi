use crate::cli::ArgParser;

use rust_htslib::{
    bam::record::Aux,
    bam::{self, IndexedReader, Read},
    bam::ext::BamRecordExtensions,
};

use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    args: ArgParser,
}

// Coverage, (hp0, hp1, hp2)
pub type DataTuple = (u64, u64, u64, [Vec<isize>; 3]);
type CoverageTuple = (u64, u64, u64); // (sub_start, sub_end, coverage)

impl BamParser {
    // Need the params to be passed here
    pub fn new(bam_name: PathBuf, reference: Option<PathBuf>, args: ArgParser) -> Self {
        let mut bam = IndexedReader::from_path(bam_name).unwrap();
        if let Some(ref_name) = reference {
            let _ = bam.set_reference(ref_name.clone());
        }
        BamParser { bam, args }
    }

    /// Extract per-base coverage/indel-delta stats for every sub-interval in
    /// a chunk, using a single BAM fetch + single pass over the reads that
    /// overlap the chunk (instead of one fetch per sub-interval).
    ///
    /// `chunk_start`/`chunk_end` should be the chunk's overall
    /// (min_start, max_end) span; `sub_intervals` the original
    /// (start, end) pairs it contains, sorted ascending by start
    /// (as produced by `chunk_regions`). They may overlap/nest.
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

        let n = sub_intervals.len();
        let mut coverage: Vec<u64> = vec![0; n];
        let mut deltas: Vec<[Vec<isize>; 3]> = (0..n)
            .map(|_| [Vec::new(), Vec::new(), Vec::new()])
            .collect();

        // Reused scratch buffers, reset per-read instead of reallocated.
        let mut diff_acc: Vec<isize> = vec![0; n];
        let mut active: Vec<usize> = Vec::with_capacity(4);

        let mut alignment = bam::Record::new();

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

            // Only sub-intervals starting strictly after ref_start can
            // possibly be spanned "from the left" by this read.
            let si_start = sub_intervals.partition_point(|&(s, _)| s <= ref_start);
            if si_start >= n {
                continue;
            }

            // Reset scratch state for this read.
            for v in diff_acc[si_start..].iter_mut() {
                *v = 0;
            }
            active.clear();
            let mut next_idx = si_start;

            let mut read_pos = ref_start;

            for &raw in alignment.raw_cigar() {
                let len = (raw >> 4) as u64;
                match raw & 0xf {
                    1 => {
                        // Ins: point event at read_pos, doesn't consume reference.
                        update_active(&mut active, &mut next_idx, sub_intervals, read_pos, n);
                        for &idx in active.iter() {
                            diff_acc[idx] += len as isize;
                        }
                    }
                    2 => {
                        // Del: point event at read_pos (before advancing).
                        update_active(&mut active, &mut next_idx, sub_intervals, read_pos, n);
                        for &idx in active.iter() {
                            diff_acc[idx] -= len as isize;
                        }
                        read_pos += len;
                    }
                    0 | 7 | 8 => {
                        read_pos += len;
                    }
                    _ => {}
                }

                if read_pos > chunk_end {
                    break; // can't reach any further sub-interval in this chunk
                }
            }

            // Commit coverage/delta for every candidate sub-interval this
            // read actually spans (read reached past its end). Ends aren't
            // guaranteed monotonic (sub-intervals may nest), so check all
            // candidates rather than breaking on the first miss.
            let mut hp: Option<usize> = None;
            for j in si_start..n {
                if sub_intervals[j].1 < read_pos {
                    let hp_val = *hp.get_or_insert_with(|| Self::extract_hp(&alignment));
                    coverage[j] += 1;
                    deltas[j][hp_val].push(diff_acc[j]);
                }
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

    /// Like `extract_reads_plup_fast`, but only counts reads spanning each
    /// sub-interval — no indel deltas, no HP splitting. One fetch + one
    /// CIGAR pass per chunk.
    pub fn extract_reads_coverage_only(
        &mut self,
        chrom: &String,
        chunk_start: u64,
        chunk_end: u64,
        sub_intervals: &[(u64, u64)],
    ) -> Vec<CoverageTuple> {
        if let Err(e) = self.bam.fetch((&chrom, chunk_start, chunk_end)) {
            panic!(
                "Unable to fetch bam {}:{}-{}\n{:?}",
                chrom, chunk_start, chunk_end, e
            )
        };

        let n = sub_intervals.len();
        let mut coverage: Vec<u64> = vec![0; n];

        let mut alignment = bam::Record::new();

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

            // Only sub-intervals starting strictly after ref_start can
            // possibly be spanned "from the left" by this read.
            let si_start = sub_intervals.partition_point(|&(s, _)| s <= ref_start);
            if si_start >= n {
                continue;
            }

            let read_pos = alignment.reference_end() as u64;

            // Commit coverage for every candidate sub-interval this read
            // spans end-to-end. Checked over all candidates (not just up
            // to the first miss) since sub-intervals may nest/overlap.
            for j in si_start..n {
                if sub_intervals[j].1 < read_pos {
                    coverage[j] += 1;
                }
            }
        }

        sub_intervals
            .iter()
            .zip(coverage)
            .map(|((s, e), cov)| (*s, *e, cov))
            .collect()
    }
}

/// Advance `next_idx`/`active` so that `active` holds exactly the indices
/// of sub-intervals whose [start, end] currently contains `pos`.
fn update_active(
    active: &mut Vec<usize>,
    next_idx: &mut usize,
    sub_intervals: &[(u64, u64)],
    pos: u64,
    n: usize,
) {
    while *next_idx < n && sub_intervals[*next_idx].0 <= pos {
        active.push(*next_idx);
        *next_idx += 1;
    }
    active.retain(|&idx| sub_intervals[idx].1 >= pos);
}
