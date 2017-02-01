`mixemt` - Deconvolving mtDNA Mixtures by Expectation Maximization
==================================================================

`mixemt` is a program for making sense of mixtures of human mitochondrial
sequences. It takes as input a reference sequence, a representation of the
mitochondrial haplogroup phylogeny from
[Phylotree.org](http://www.phylotree.org/), and a BAM file containing reads
mapped to the reference sequence. `mixemt` scans each mapped read for variants
described in Phylotree and estimates the probability of the read originating
from each candidate haplotype (haplogroup). This information is used to
estimate the mixture proportions for the sample and to identify the most likely
haplogroups contributing to the sample. The program produces as output a table,
written to standard output, that describes the detected haplogroups and their
estimated percent contribution to the sample. The program can also produce new
BAM files of input fragments partitioned by the most likely contributor of
origin, tab-delimited files containing statistics for each reference position,
and verbose output detailing the results from each step.

## Requirements

`mixemt` is written in Python and requires a few additional packages:

* Numpy
* Biopython
* pysam

R and ggplot2 are required to use the plotting scripts included in the
directory `plot/`.

## Usage

```
mixemt [options] <ref_seq.fa> <phylotree.csv> <reads.bam>
```

## Preparing input sequences

`mixemt` takes as input a reference sequence (`ref/`), a representation of
the phylogeny from Phylotree (`phylotree/`), and a BAM file containing
alignments of sequences to the reference. The Reconstructed Sapiens Reference
Sequence (RSRS) is used the reference sequence, although reads may be mapped to
the Revised Cambridge Reference Sequence (rCRS) as both share a common
coordinate system.

Aligned sequences are provided to `mixemt` as a BAM file. Alignments should be
filtered prior to input using `samtools` or a similar program.  Internally,
`mixemt` treats alignments sharing the same query template name as a single
fragment and does not check read pairs for orientation, insert size, unaligned
segments, or secondary alignments. The options `-q` and `-Q` allow for basic
filtering of alignments and base calls based on map quality and base quality
scores.

Since the mitochondrial genome is circular, it may be advantageous to align
reads across the junction of the start and end of the reference sequence,
split the reads, and correct the mapping coordinates. A basic implementation of
these steps can be found in the repository
[circ\_aln](https://github.com/svohr/circ_aln).


## Options

