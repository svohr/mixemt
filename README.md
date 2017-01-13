`mixemt` - Deconvoluting mtDNA Mixtures by Expectation Maximization
==================================================================

`mixemt` is a program for making sense of mixtures of human mitochondrial
sequences. It takes as input a reference sequence (usually the Reconstructed
Sapiens Reference Sequence), a representation of the mitochondrial haplogroup
phylogeny from [Phylotree.org], and a BAM file containing reads mapped to the
reference sequence. `mixemt` scans each mapped read for variants described in
Phylotree and estimates the probablity of the read originating from every
haplogroup.  This information is used to estimate the mixture proportions for
the sample and to identify the most likely haplogroups contributing to the
sample. The program produces as output a table, written to standard output,
that describes the detected haplogroups and their percent contribution to the
sample. The program can also produce new BAM files of input fragments
partitioned by the most likely contributor of origin, tab-delimited files
containing statistics for each reference position, and verbose output detailing
the results from each step.
