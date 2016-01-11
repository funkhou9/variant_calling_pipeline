# Scott Funkhouser's Variant Calling Pipeline

### Used to isolate variants from high-throughput sequencing data

> This pipeline has been applied to the discovery of candidate RNA editing sites. In this
> case, resulting .vcf files are processed with [`editTools`](https://github.com/funkhou9/editTools)

## Table of Contents

1. [Read Trimming](#read-trimming)
2. [Mapping](#mapping)
3. [Merging and Indexing](#merging-and-indexing)
4. [Variant Calling](#variant-calling)

## Read Trimming

Raw sequencing reads are trimmed with [`condetri`](https://github.com/linneas/condetri), which
allows quality filtering at the 3' ends of paired or single-end sequencing data.

### Example condetri usage with paired-end data:

```sh
perl condetri_v2.2.pl -fastq1=<lane1_read1.fastq> -fastq2=<lane1_read2.fastq> -sc=33 -minlen=75 -prefix=<lane1> -pb=fq
```

`-minlen` sets the minimum length of reads after trimming.  
`-sc` is used to provide the scoring scheme of the base qualities. Modern Illumina data uses 33.  
`-prefix` sets the prefix of the output files.  
`-pb=fq` prints low quality reads that were removed during the trimming process to a fastq
file for diagnostic purposes.  
`<lane1_read1.fastq>` and `<lane1_read2.fastq>` represent paired-end reads from the same lane.

### Output

```
 lane1_trim1.fastq             File(s) with the trimmed reads (one file for
 lane1_trim2.fastq              single-end data, three for paired-end data, where
 lane1_trim_unpaired.fastq      the last file includes reads from the two input
                                files whose read pair had too poor quality.
 lane1.stats                   Includes basic statistics in columns.
```

In the case of paired-end data from multiple lanes, subsequent lanes are processed the same as the first.

> Manual trimming of the 5' ends of reads is also possible with `condetri`. This is used to process RNA
> sequencing data for the purposes of RNA editing detection. In this case, the additional parameter `-cutfirst=6`
> is used to trim the first 6 bases from the 5' ends of reads. An explanation for doing so can be found in
> [Bass et al. 2012](http://www.nature.com/nbt/journal/v30/n12/full/nbt.2452.html).

## Mapping

For each lane, trimmed paired and unpaired reads in files `_trim1.fastq`, `_trim2.fastq`, and `_trim_unpaired.fastq`
are mapped to the reference genome. Depending on the source material (DNA or RNA), a different mapper may be used.

### DNA (Whole genome sequence)

Trimmed WGS reads are mapped using [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

### Example bowtie2 v usage with paired and unpaired WGS data:

```sh
bowtie2 -x <reference> -p 7 -X 1000 -1 <_trim1.fastq> -2 <_trim2.fastq> -U <_trim_unpaired.fastq> -S <out.sam>
```

`-x` supplies the basename of the index for the reference genome  
`-p` sets the number of processors to be used in parallel (assuming multi-processor machine)  
`-X` sets the maximum fragment length for valid paired-end alignments  
`-1`, `-2`, and `-U` supply the first pair, second pair and unpared reads, respectively

### Output

```
out.sam		"Sequence alignment file" (SAM).
			Contains alignment information for each read.
```

For more information on SAM format, click [here](http://samtools.github.io/hts-specs/SAMv1.pdf).

### RNA (RNAseq)

Alternatively, for sequencing reads generated from cDNA, the splice-site aware aligner [`TopHat`](https://ccb.jhu.edu/software/tophat/index.shtml) is used.

### Example TopHat v usage with paired and unpaired RNASeq data:

```sh

```


