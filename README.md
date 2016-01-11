# Scott Funkhouser's Variant Calling Pipeline

### Used to isolate variants from high-throughput sequencing data

> This pipeline has been applied to the discovery of candidate RNA editing sites. In this
> case, resulting VCF files are processed with [`editTools`](https://github.com/funkhou9/editTools). Optional notes such as these pertaining to RNA editing can be found within blockquotes throughout this document.

## Table of Contents

1. [Read Trimming](#read-trimming)
2. [Mapping](#mapping)
3. [Converting, Merging, Sorting and Indexing](#converting,-merging,-sorting,-and-indexing)
4. [Variant Calling](#variant-calling)

## Read Trimming

Raw sequencing reads are trimmed with [`condetri`](https://github.com/linneas/condetri), which
allows quality filtering at the 3' ends of paired or single-end sequencing data.

#### Example condetri usage with paired-end data:

```sh
perl condetri_v2.2.pl -fastq1=<lane1_read1.fastq> -fastq2=<lane1_read2.fastq> -sc=33 -minlen=75 -prefix=<lane1> -pb=fq
```

`-minlen` sets the minimum length of reads after trimming.  
`-sc` is used to provide the scoring scheme of the base qualities. Modern Illumina data uses 33.  
`-prefix` sets the prefix of the output files.  
`-pb=fq` prints low quality reads that were removed during the trimming process to a fastq
file for diagnostic purposes.  
`<lane1_read1.fastq>` and `<lane1_read2.fastq>` represent paired-end reads from the same lane.

#### Output

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

#### DNA (Whole genome sequence)

Trimmed WGS reads are mapped using [`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

#### Example Bowtie 2 v2.2.1 usage with paired and unpaired WGS data:

```sh
bowtie2 -x <reference> -p 7 -X 1000 -1 <_trim1.fastq> -2 <_trim2.fastq> -U <_trim_unpaired.fastq> -S <out.sam>
```

`-x` supplies the basename of the index for the reference genome  
`-p` sets the number of processors to be used in parallel (assuming multi-processor machine)  
`-X` sets the maximum fragment length for valid paired-end alignments  
`-1`, `-2`, and `-U` supply the first pair, second pair and unpared reads, respectively

#### Output

```
out.sam		"Sequence alignment file" (SAM).
			Contains alignment information for each read.
```

For more information on SAM format, click [here](http://samtools.github.io/hts-specs/SAMv1.pdf).

#### RNA (RNAseq)

Alternatively, for sequencing reads generated from cDNA, the splice-site aware aligner [`TopHat`](https://ccb.jhu.edu/software/tophat/index.shtml) is used.

#### Example TopHat v2.0.12 usage with paired and unpaired strand-specific RNASeq data:

```sh
tophat -o <out> -p 7 --mate-inner-dist 400 --mate-std-dev 100 --library-type "fr-firststrand" <reference>  <_trim1.fastq>,<_trim_unpaired.fastq> <_trim2.fastq>
```

`-o` supplies the basename of the SAM output  
`-p` sets the number of processors to be used in parallel (assuming multi-processor machine)  
`--mate-inner-dist` sets the expected inner distance between mate pairs  
`--library-type` sets the protocol used to generate strand-specific libraries. `"fr-firststrand"` is the protocol used by MSU's genomics core

#### Output

```
out.bam		"Binary alignment file" (BAM)
```
## Converting, Merging, Sorting and Indexing

Prior to using alignment files for variant calling, files may need converting, libraries need merging, and merged alignments need sorting and indexing. These steps are acheived through a combination of software, [`Samtools`](http://www.htslib.org) and [`Picard`](hhttps://github.com/broadinstitute/picard).

#### DNA (Whole genome sequence)

For each library, the output of Bowtie 2 (SAM format) is converted to BAM format using Samtools with:

```sh
samtools view -bS <out.sam> > <out.bam>
```
`-b` sets the output as BAM format  
`-S` specifies the input as SAM format

Libraries are then merged with:

```sh
samtools merge <out_merged.bam> <out1.bam> <out2.bam> ...
```

Merged libraries are then sorted using:

```py
samtools sort <out_merged.bam> <out_merged_sorted>
```

The merged and sorted output is finally indexed using:

```sh
samtools index <out_merged_sorted.bam>
```

> Variant calling may require conservative measures to filter out reads that do not "uniquely map". For the purpose of RNA editing detection, "uniquely mapped" reads are extracted from the merged, sorted BAM file by grepping reads without the Bowtie 2 tag `XS:i:`.

#### RNA (RNAseq)

The output from TopHat will already be in BAM format, but requires merging libraries using the same `samtools merge` command above.

> In some instances, such as simultaneous RNA and DNA variant calling for the purposes of RNA editing detection, the merged RNA alignment files may need to be reordered to match the contig order of both the DNA alignment files and reference file. After building [`Picard`](https://github.com/broadinstitute/picard), this can be done using 
> ```
> java -jar $PICARD/ReorderSam.jar R=<reference.fa> I=<out_merged.bam> O=<out_merged_reordered.bam>
> ```

> If "uniquely mapped" reads are desired, uniquely mapped RNASeq reads can be obtained by grepping reads that contain the `NH:i:1` tag.

> In order to separate RNASeq reads coming from plus-strand transcripts vs minus-strand transcripts (again important for RNA editing detection) reads containing the `XS:A:+` tag must be separated from those with the `XS:A:-` tag.

Indexing of merged files is accomplished using the same `samtools index` command as above.

## Variant Calling

Calling variants is done using the algorithm detailed in [Heng Li 2011](http://bioinformatics.oxfordjournals.org/content/early/2011/09/08/bioinformatics.btr509.abstract), implemented in [`Samtools`](http://www.htslib.org/doc/samtools.html) and [`bcftools`](http://www.htslib.org/doc/bcftools.html).

`samtools mpileup` is first used to calculate genotype likelihoods for each genomic position, then `bcftools call` is used to obtain posterior genotype calls in variant call format (VCF).

#### Example of variant calling with either RNAseq or WGS data

```sh
samtools mpileup -f <reference.fa> -C50 -E -Q25 -ug -t DP,DV <out_merged_sorted.bam> | bcftools call -O v -m -v > out.vcf
```
##### samtools mpileup v1.0

`-f` supplies reference file name  
`-C` sets coefficient for downgrading [mapping qualities](http://maq.sourceforge.net/qual.shtml) (MQ) for reads containing excessive mismatches. Recommended value is 50 for BWA.  
`-E` sets computation of base alignement qualities and ignores standard base qualities. See [Heng Li 2011](http://www.ncbi.nlm.nih.gov/pubmed/21320865)  
`-Q` sets minimum base quality or base alignment quality for a base to be considered
`-ug` sets format of output to uncompressed bcf (suggested for piping)  
`-t` used to specify desired fields in output


##### bcftools call v1.0

`-O` sets output type (v for VCF)  
`-m` sets the 'multiallelic caller', an alternative model that overcomes limitations of previous models  
`-v` sets output to only contain variant sites rather than statistics for all genomic positions

> For RNA editing detection, both DNA and RNA alignments are used simultanously in the variant calling step, specifying both files with:
> 
> ```sh
> samtools mpileup -f <reference.fa> -C50 -E -Q25 -ug -t DP,DV <out_merged_sorted_dna.bam> <out_merged_reordered_rna.bam> | bcftools call -O v -m -v > out.vcf
> ```
> out.vcf is then processed with [`editTools`](https://github.com/funkhou9/editTools).

 

 







