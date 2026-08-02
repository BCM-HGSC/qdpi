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
