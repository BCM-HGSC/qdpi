Quick Delta and PIleups
=======================

Over a set of bed regions, extract the coverage of reads spanning the region and their length delta by HP tag.

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
of allele length over a region. 

In order to ensure that changes observed over a bed-region actually belong into said bed-region, alignments over
flanking regions around the bed-region are also checked. First, only alignments spanning `--flank` base-pairs 
upstream/downstream of the bed-region are considered. Second, only alignments which have matching bases to at least
`--min-anchor` percent of the bed-region's flank are considered. The defaults of `--flank 200` and `--min-anchor 0.90`
roughly enforces a 90% similarity between read and flanking-regions. Note that this similarity is an overestimate since
insertions are not counted against the anchor sequence. Furthermore, when two bed-regions are within `--flank` and one
has a larger deletion, this will prevent the other from collecting coverage.


