use serde::Serialize;
use std::io::{BufWriter, Write};
use std::fs::File;

fn mean_std(data: &[isize]) -> (f64, f64) {
    let n = data.len() as f64;
    let sum: isize = data.iter().sum();
    let mean = sum as f64 / n;
    let variance = data.iter().map(|&x| ((x as f64) - mean).powi(2)).sum::<f64>() / n;
    let std_deviation = variance.sqrt();
    (mean, std_deviation)
}

#[derive(Debug, Serialize)]
pub struct Locus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: usize,
    pub n_clusters: usize,
    pub cluster_score: f32,
    pub clusters: Vec<Vec<usize>>,
    #[serde(skip_serializing)]
    pub kfeats: Vec<Vec<f32>>,
    #[serde(skip_serializing)]
    pub read_lengths: Vec<usize>,
    pub read_seqs: Vec<Vec<u8>>,
    //pub read_lengths: HashMap<usize, usize>,
}

// chrom    start   end coverage    n_clusters  list_of_vaf list_of_deltas
impl Locus {
    pub fn new(chrom: &String, start: u64, end: u64) -> Self {
        Locus {
            chrom: chrom.to_string(),
            start,
            end,
            n_reads: 0,
            n_clusters: 0,
            cluster_score: 0.0,
            read_lengths: vec![],
            clusters: vec![],
            kfeats: vec![],
            read_seqs: vec![],
        }
    }
    
    pub fn make_verbose_data(&self, file: &mut BufWriter<File>) -> std::io::Result<()> {
        let reg_key = format!("{}:{}-{}", self.chrom, self.start, self.end);
        for (idx, m_clust) in self.clusters.iter().enumerate() {
            for read_index in m_clust.iter().filter(|&c| !self.read_seqs[*c].is_empty()) {
                write!(file, ">{}_{}\n{}\n",
                reg_key, idx, String::from_utf8(self.read_seqs[*read_index].clone()).expect("Our bytes should be valid utf8"))?;

            }
        }
        Ok(())
    }

    pub fn make_data(&self, file: &mut BufWriter<File>) -> std::io::Result<()> {
        write!(file, "{}\t{}\t{}\t{}\t{}\t{}", self.chrom, self.start, self.end, self.n_reads, self.n_clusters, self.cluster_score)?;
        
        let span: isize = (self.end as isize) - (self.start as isize);
        for m_clust in &self.clusters {
            let cluster_read_lengths = m_clust.iter().map(|&read_index| self.read_lengths[read_index] as isize - span).collect::<Vec<isize>>();
            let (mean_length, std_deviation) = mean_std(&cluster_read_lengths);
            let count = cluster_read_lengths.len();

            write!(file,
                    "\t{:.2},{:.2},{}",
                    mean_length,
                    std_deviation,
                    count
            )?;
        }
        write!(file, "\n")
    }
}
