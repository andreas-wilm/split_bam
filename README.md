Split BAM
=========

Splits your BAM/SAM stream into one file per region given in bed file
(plus one for unaligned reads).

Similar to [bamUtils' splitBam](http://genome.sph.umich.edu/wiki/BamUtil:_splitBam) but no read group
dependency

Main purpose is to take aligner output and split it into different
files, which can be processed in parallel by tools that work on the
same region (e.g. reference).

Replaces my own [split_bam_by_chr](https://github.com/andreas-wilm/split_bam_by_chr)

# Compilation

- You will have to set HTSLIBDIR. 
- Default target is 'static'. The other valid target is 'dynamic'
- Example: `HTSLIBDIR=/opt/local/htslib-1.3/ make dynamic`
