Quick Delta and PIleups
=======================

Over a set of bed regions, extract the coverage, read-length delta, and pileup variant positions.

Output is a tsv bed-like with columns:
* chromosome
* start
* end
* data tuple with three values
  * coverage
  * list of [position, variant]
  * list of read-length deltas

The position is relative to start (e.g. `start + position = genomic position`). Variant is a negative number for
deletions or inserted sequences. SNPs are not recorded.
