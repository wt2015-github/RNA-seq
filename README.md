# RNA Sequencing (RNA-seq) Gene Expression Data Extraction

An automated Unix shell pipeline for quality control (QC) and extracting gene-expression values from RNA sequencing (RNA-seq) raw data files (*.fastq* or *fastq.gz*)

Original script author: Dr. Ting Wang (@wt2015-github via [repo](github.com/wt2015-github/RNA-seq)); Adapted and re-tested by: David Chen (@ydavidchen via [public fork](github.com/ydavidchen/RNA-seq))

## Overview

### Credit

Please refer to the original author's [repository](github.com/wt2015-github/RNA-seq). This [remote GitHub branch](github.com/ydavidchen/RNA-seq) has been adapted to and optimized for extracting gene expression data.

Please also refer to the relevant citations for open-source software used:
* [RSEM](https://www.ncbi.nlm.nih.gov/pubmed/21816040)
* [samtools](www.ncbi.nlm.nih.gov/pubmed/19505943)
* [STAR](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
* [multiqc](www.ncbi.nlm.nih.gov/pubmed/27312411)
* [fastqc (no published citation)](www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [trim_galore (no published citation)](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

### Scripts and Files

* *Dependencies.sh* includes commands and instructions for installing relevant open-source bioinformatics software, and setting executable aliases. This script is to be executed in real-time *interactively*. Adaptation to user's own computer machine is most likely necessary.
* *run_RNAseq_STAR.RSEM.sh* includes QC, trimming, contamination screening, mapping, gene/isoform expression quantification, and summary steps (adapted from the [original author](github.com/wt2015-github/RNA-seq)). This script is to be executed as its entirety in command line.

Upon successful completion of this Unix shell pipeline, you may use gene expression matrices exported for downstream analyses, which are often performed in the [R statistical programming environment](www.r-project.org).

# Usage

## Syntax

`run_RNAseq_STAR.RSEM.sh -g hg19 <-p> -i <path_of_inputs> -o <path_of_outputs> -t <threads>`

## Arguments

* *-g*: Set reference genome. Defaults to *hg19*.
* *-p*: Set this for paired-end sequencing design. Defaults to single-end design.
* *-i*: Directory of input *fastq* or *fastq.gz* files.
* *-o*: Directory for output files.
* *-t*: Set average number of threads for **each** sample in parallel analysis. Must be an integer value. Defaults to $1$.

## Configurations

The paths of some genome files and tools need to be modified accordingly:
* Download, and then modify, the paths of *genome STAR index*, *ref*, *fa*, *gtf* and *refSeq*.
* Make sure `fastqc`, `multiqc`, `trim_galore`, `fastq_screen`, `STAR`, `RSEM`, and `samtools` are installed and can be executed by just typing their names. Otherwise, you need to modify the paths of these tools in*run_RNAseq_STAR.RSEM.sh*.
* Install `Picard` tools and set an alias (environment variable) named "`PICARD`" for the path of *picard.jar* script. Otherwise, you can manually change "`$PICARD`" in this pipeline script to the full path of *picard.jar* script.
* Add *summarize_RNAseq_STAR.RSEM.pl* to environment `PATH` or put it in */usr/local/bin/*. Otherwise, modify the path of this Perl script in this *run_RNAseq_STAR.RSEM.sh* pipeline.

## Outputs

* *qc* directory: QC results (e.g. *multiqc* summary html report)
* *process* directory consisting of:
  - Summary table
  - Gene or isoform expression count/FPKM/TPM matrix
  - Sub-directories of analysis results for each sample

# Example Data

## Data Set

A data set can be found here. Please download the paired-end SRR runs and store in the user-specified directory. As a reminder, storing as the compressed *.fastq.gz* format is fine.

[Paired-end Illumina HiSeq2000 breast cancer MCF7 cell line](https://www.ebi.ac.uk/ena/data/view/SRR1021668)

## Procedure

* Download all dependencies by interactively running *Dependency *
* Download human hg19 genome, gene, and transcriptome annotations via [Ensembl ftp site](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/) and [UCSC Table Browser](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables)

# Coming soon

Example optimized for remote computation on Linux cluster
