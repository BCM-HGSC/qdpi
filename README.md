Quick Delta and PIleups
=======================

Over a set of bed regions, extract the coverage, read-length delta, and pileup variant positions.

Output is a tsv bed-like with columns:
* chromosome
* start
* end
* coverage
* hp0 - comma-separated list of haplotype 0 read lengths.
* hp1 - comma-separated list of haplotype 1 read lengths.
* hp2 - comma-separated list of haplotype 2 read lengths.

Nulls (e.g. no unphased `hp0` reads) will be `.`

Details
=======

Aligned read CIGAR fields from a BAM/CRAM are parsed. Subsequences spanning the bed coordinates are selected and every
insertion/deletion operation summed to create a "Read Delta". This read delta represents the relative increase/decrease
of an allele over a region. 

In order to ensure that changes observed over a bed-region actually belong into said bed-region, alignments over
flanking regions around the bed-region are also checked. First, only alignments spanning `--anchor` base-pairs 
upstream/downstream of the bed-region are considered. Second, only alignments with fewer then `--max-edits` edit 
operations over the bed-region are considered. The defaults of `--anchor 200` and `--max-edits 10` roughly 
enforces a maximum flanking-region error rate of 5%
