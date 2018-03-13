# RNA-seq analysis
A shell pipeline for automatically analyzing all RNA-seq data in a folder based on STAR-RSEM and GATK

## Introduction

*run_RNAseq_STAR.RSEM.sh* includes QC, trimming, contamination screening, mapping, gene/isoform expression quantification, and summary steps.

*callSNP_RNAseq.sh* is for SNP calling using the RNA-seq analysis results for each sample.

Downstream analyses, such as different expression analysis and pathway enrichment analaysis, are case specific, so are not included in this pipeline.

## Configuration
The paths of some genome files and tools need to be modified accordingly:
* Modify the paths of genome STAR index and ref, fa, gtf, refSeq, and vcf (if use SNP calling) files.
* Make sure **fastqc**, **multiqc**, **trim_galore**, **fastq_screen**, **STAR**, **RSEM**, **samtools** are installed and can be ran by just typing their names, otherwise modify the paths of these tools in this pipeline script.
* Install **Picard** tools and set an environment variable/alias named "PICARD" for the path of picard.jar script, otherwise manually change "**$PICARD**" in this pipeline script to the full path of picard.jar script. The same with **GATK**, if use the SNP calling script. 
* Add summarize_RNAseq_STAR.RSEM.pl to environment PATH or put it in /usr/local/bin/, otherwise modify the path of this perl file in this pipeline script.

## Usage
```
run_RNAseq_STAR.RSEM.sh g <hg19, mm10 or rn6> <-p> -i <path_of_inputs> -o <path_of_outputs> -t <threads>

callSNP_RNAseq.sh -g <hg19> -i <input_folder> -t <threads>
```

## Arguments
* *-g*: set hg19, mm10 or rn6 as reference genome, default is hg19.
* *-p*: set this for paired-end data, default is for single-end data
* *-i*: directory of input fastq or fastq.gz files.
* *-o*: directory for output files.
* *-t*: set average number of threads for **each** sample in parallel analysis, must be integer, default is 1.

## Outputs
*run_RNAseq_STAR.RSEM.sh* generates:
* A **qc** folder including QC results.
* A **process** folder including a summary table, gene/isoform expression count/FPKM/TPM tables, and subfolders of analysis results for each sample.
*callSNP_RNAseq.sh* generates SNP calling vcf files in each sample folder.

## Contact
[Ting Wang](http://wt2015-github.github.io/) ([email](wang9ting@gmail.com)).
