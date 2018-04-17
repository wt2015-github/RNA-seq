#!/bin/bash
#Author: T. Wang, 2017-03

manual="
Shell script for automatically analyzing all RNA-seq data in a folder based on STAR-RSEM.
QC with  at the beginning and use multiqc to summary qc figures.
Note RSEM --star can not implement STAR for all samples with a shared memory of genome, so separate the two steps!
Usage:
    run_RNAseq_RSEM.sh -g <hg19> <-p> -i <path_of_inputs> -o <path_of_outputs> -t <threads>
Options:
    -g    set hg19 as reference genome, default is hg19
    -p    set this for paired-end data, default is for single-end data
    -i    directory of input fastq or fastq.gz files
    -o    directory for output files, please use full pathway
    -t    average number of threads for each sample, must be integer, default is 1
"

if [[ $# -le 0 ]]; then
    echo "${manual}" >&2
    exit 1
fi

while getopts :g:pi:o:t: ARGS
do
case $ARGS in
    g)
        genome=$OPTARG
        ;;
    p)
        datatype="PE"
        ;;
    i)
        pathin=$OPTARG
        ;;
    o)
        pathout=$OPTARG
        ;;
    t)
        threads=$OPTARG
        ;;
    :)
        echo "no value for option: $OPTARG"
        echo "${manual}" >&2
        exit 1
        ;;
    *)
        echo "unknow option: $OPTARG"
        echo "${manual}" >&2
        exit 1
        ;;
esac
done

cd .
echo "Current working directory is:"
echo $PWD

echo ${genome:="hg19"} >/dev/null
echo ${datatype:="SE"} >/dev/null
echo ${threads:=1} >/dev/null
echo ${pathin:="nopathin"} >/dev/null
echo ${pathout:="nopathout"} >/dev/null

if [[ ${pathin} == "nopathin" ]]; then
    echo "=== Please set directory of inputs === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathin} ]]; then
    echo "=== Input directory does not exist! === ${manual}" >&2
    exit 1
fi

if [[ ${pathout} == "nopathout" ]]; then
    echo "=== Please set directory of outputs === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathout} ]]; then
    echo "=== Output directory does not exist! === ${manual}" >&2
    exit 1
fi

(( ${threads} )) 2>/dev/null
if [[ $? != 0 || ! ${threads} -ge 1 ]]; then
    echo "=== Thread number must be a positive integer! === ${manual}" >&2
    exit 1
fi
echo ""
echo "=== Set ${threads} threads for each sample in parallel analysis ==="

if [[ ${genome} == "hg19" ]]; then
    echo "=== set reference genome: hg19 ==="
    genome_RSEM_STAR_index="/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/RSEM_STAR_Index"
    genome_RSEM_STAR_ref="/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/RSEM_STAR_Index/human_ref"
    genome_fa="/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    genome_gtf="/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
    genome_refSeq_Bed="/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_RefSeq.bed"
else
    echo "=== Supported reference genomes is hg19 === ${manual}" >&2
    exit 1
fi

fastq_files=(`find ${pathin} -maxdepth 1 -name "*.fastq*" -type f | sort`)
if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "=== No .fastq or .fastq.gz file found in input directory! === ${manual}" >&2
    exit 1
else
    echo "=== There are ${#fastq_files[@]} raw fastq files in input directory ==="
fi

if [[ ${datatype} == "PE" ]]; then
    echo "=== Quality control with FastQC ==="
    mkdir ${pathout}/fastqc
    for file in ${fastq_files[@]}; do fastqc $file -t ${threads} -o ${pathout}/fastqc/ & done
    wait
    fastqc_files=(`ls ${pathout}/fastqc/*fastqc.zip`)
    if [[ ${#fastq_files[@]} -ne ${#fastqc_files[@]} ]]; then
        echo "=== Some samples did not pass fastqc! ===" >&2
        exit 1
    fi
    multiqc -o ${pathout}/qc/ ${pathout}/fastqc/
    mv ${pathout}/fastqc ${pathout}/qc/
    echo "=== Load-genome into shared memory for parallel alignment ==="
    mkdir ${pathout}/tmp
    STAR --genomeDir ${genome_RSEM_STAR_index} --genomeLoad LoadAndExit --outFileNamePrefix ${pathout}/tmp/
    rm -rdf ${pathout}/tmp
    echo "=== Parallel-analyzing paired-end RNA-seq data, see log.all.txt in each subfolder ==="
    mkdir ${pathout}/process
    for ((i=0; i<${#fastq_files[@]}; i+=2)); do
        subfolder=$(echo ${fastq_files[$i]} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)
        mkdir ${pathout}/process/${subfolder}
        echo "=== analyzing for sample ${subfolder} ==="
        {
            echo "=== trimming data and QC for sample ${subfolder} ==="
            trim_galore ${fastq_files[$i]} ${fastq_files[$i+1]} --length 25 --suppress_warn --paired --gzip --o ${pathout}/process/${subfolder}/
            echo "=== screening contamination for sample ${subfolder} ==="
            fastq_screen ${pathout}/process/${subfolder}/*.fq.gz --subset 1000000 --threads ${threads} --quiet
            echo "=== mapping to reference genome with STAR for sample ${subfolder} ==="
            STAR --genomeDir ${genome_RSEM_STAR_index} --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --runThreadN ${threads}  --genomeLoad LoadAndKeep  --outSAMtype BAM Unsorted  --quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted  --outFileNamePrefix ${pathout}/process/${subfolder}/  --readFilesCommand zcat  --readFilesIn ${pathout}/process/${subfolder}/*.fq.gz
            echo "=== assembling transcripts with RSEM based on STAR alignments for sample ${subfolder} ==="
            rsem-calculate-expression --strandedness reverse -p ${threads} -q --time --paired-end --alignments ${pathout}/process/${subfolder}/Aligned.toTranscriptome.out.bam ${genome_RSEM_STAR_ref} ${pathout}/process/${subfolder}/rsem
            echo "=== sort genome alignments by coordinate for sample ${subfolder} ==="
            samtools sort -@ ${threads} -m 1G -o ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam ${pathout}/process/${subfolder}/Aligned.out.bam
            samtools index ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam
            echo "=== annotating bam for sample ${subfolder} ==="
            read_distribution.py -i ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam -r ${genome_refSeq_Bed} > ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam.distribution
            mkdir ${pathout}/process/${subfolder}/tmp_java
            echo "=== marking duplicates for sample ${subfolder} ==="
            java -Xmx10g -XX:ParallelGCThreads=${threads} -jar $PICARD MarkDuplicates I=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam O=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam M=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=${pathout}/process/${subfolder}/tmp_java/
            echo "=== getting flag states for sample ${subfolder} ==="
            samtools flagstat ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam > ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.flagstat
            rm -rdf ${pathout}/process/${subfolder}/tmp_java
            rm ${pathout}/process/${subfolder}/Aligned.out.bam
        } > ${pathout}/process/${subfolder}/log.all.txt 2>&1 &
    done
wait
elif [[ ${datatype} == "SE" ]]; then
    echo "=== Quality control with FastQC ==="
    mkdir ${pathout}/fastqc
    for file in ${fastq_files[@]}; do fastqc $file -t ${threads} -o ${pathout}/fastqc/ & done
    wait
    fastqc_files=(`ls ${pathout}/fastqc/*fastqc.zip`)
    if [[ ${#fastq_files[@]} -ne ${#fastqc_files[@]} ]]; then
        echo "=== some samples not pass fastqc!!===" >&2
        exit 1
    fi
    multiqc -o ${pathout}/qc/ ${pathout}/fastqc/
    mv ${pathout}/fastqc ${pathout}/qc/
    echo "=== load genome into shared memory for parallel alignment ==="
    mkdir ${pathout}/tmp
    STAR --genomeDir ${genome_RSEM_STAR_index} --genomeLoad LoadAndExit --outFileNamePrefix ${pathout}/tmp/
    rm -rdf ${pathout}/tmp
    echo "=== parallel analyzing single-end RNA-seq data, see log.all.txt in each subfolder ==="
    mkdir ${pathout}/process
    for ((i=0; i<${#fastq_files[@]}; i++)); do
        subfolder=$(echo ${fastq_files[$i]} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)
        mkdir ${pathout}/process/${subfolder}
        echo "=== analyzing for sample ${subfolder} ==="
        {
            echo "=== trimming data and QC for sample ${subfolder} ==="
            trim_galore ${fastq_files[$i]} --length 25 --suppress_warn --gzip --o ${pathout}/process/${subfolder}/
            echo "=== screening contamination for sample ${subfolder} ==="
            fastq_screen ${pathout}/process/${subfolder}/*.fq.gz --subset 1000000 --threads ${threads} --quiet
            echo "=== mapping to reference genome with STAR for sample ${subfolder} ==="
            STAR --genomeDir ${genome_RSEM_STAR_index} --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --runThreadN ${threads}  --genomeLoad LoadAndKeep  --outSAMtype BAM Unsorted  --quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted  --outFileNamePrefix ${pathout}/process/${subfolder}/  --readFilesCommand zcat  --readFilesIn ${pathout}/process/${subfolder}/*.fq.gz
            echo "=== assembling transcripts with RSEM based on STAR alignments for sample ${subfolder} ==="
            rsem-calculate-expression --strandedness reverse -p ${threads} -q --time --alignments ${pathout}/process/${subfolder}/Aligned.toTranscriptome.out.bam ${genome_RSEM_STAR_ref} ${pathout}/process/${subfolder}/rsem
            echo "=== sort genome alignments by coordinate for sample ${subfolder} ==="
            samtools sort -@ ${threads} -m 1G -o ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam ${pathout}/process/${subfolder}/Aligned.out.bam
            samtools index ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam
            echo "=== annotating bam for sample ${subfolder} ==="
            read_distribution.py -i ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam -r ${genome_refSeq_Bed} > ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam.distribution
            mkdir ${pathout}/process/${subfolder}/tmp_java
            echo "=== marking duplicates for sample ${subfolder} ==="
            java -Xmx10g -XX:ParallelGCThreads=${threads} -jar $PICARD MarkDuplicates I=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam O=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam M=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=${pathout}/process/${subfolder}/tmp_java/
            echo "=== getting flag states for sample ${subfolder} ==="
            samtools flagstat ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam > ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.flagstat
            rm -rdf ${pathout}/process/${subfolder}/tmp_java
            rm ${pathout}/process/${subfolder}/Aligned.out.bam
        } > ${pathout}/process/${subfolder}/log.all.txt 2>&1 &
    done
wait
else
    echo "=== Wrong input of data type (PE or SE) === ${manual}" >&2
    exit 1
fi

echo "=== Analysis finished! Releasing genome from shared memory... ==="
mkdir ${pathout}/tmp
STAR --genomeDir ${genome_RSEM_STAR_index} --genomeLoad Remove --outFileNamePrefix ${pathout}/tmp/
rm -rdf ${pathout}/tmp

# Generate count table
echo "=== Generating count table for all samples with featureCounts... ==="
if [[ ${datatype} == "PE" ]]; then
    featureCounts -p -C -s 2 -T 10 -Q 10 -a ${genome_gtf} -o ${pathout}/process/counts.txt ${pathout}/process/*/Aligned.sortedByCoord.out.bam
else
    featureCounts -s 2 -T 10 -Q 10 -a ${genome_gtf} -o ${pathout}/process/counts.txt ${pathout}/process/*/Aligned.sortedByCoord.out.bam
fi
sed -i '1d' ${pathout}/process/counts.txt
cut -f 1,7- ${pathout}/process/counts.txt > ${pathout}/process/counts2.txt
mv ${pathout}/process/counts2.txt ${pathout}/process/counts.txt

# Generate expected count table
echo "=== Generating expected count table... ==="
rsem-generate-data-matrix ${pathout}/process/*/rsem.genes.results > ${pathout}/process/rsem_expected_count_genes.txt
rsem-generate-data-matrix ${pathout}/process/*/rsem.isoforms.results > ${pathout}/process/rsem_expected_count_isoforms.txt

# Generate FPKM table
echo "=== Generating FPKM table for all samples... ==="
folders=(`ls -d ${pathout}/process/*/ | rev | cut -d '/' -f2 | rev`)
mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,7 ${pathout}/process/${folders[$i]}/rsem.genes.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/fpkm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/fpkm_all.txt > ${pathout}/process/tmp/fpkm_table.txt
headarr=('gene_id' 'transcript_id(s)' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/fpkm_table.txt > ${pathout}/process/rsem_fpkm_genes.txt
rm -rdf ${pathout}/process/tmp

mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,7 ${pathout}/process/${folders[$i]}/rsem.isoforms.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/fpkm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/fpkm_all.txt > ${pathout}/process/tmp/fpkm_table.txt
headarr=('transcript_id' 'gene_id' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/fpkm_table.txt > ${pathout}/process/rsem_fpkm_isoforms.txt
rm -rdf ${pathout}/process/tmp

# Generate TPM table
echo "=== Generating TPM table for all samples... ==="
folders=(`ls -d ${pathout}/process/*/ | rev | cut -d '/' -f2 | rev`)
mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,6 ${pathout}/process/${folders[$i]}/rsem.genes.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/tpm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/tpm_all.txt > ${pathout}/process/tmp/tpm_table.txt
headarr=('gene_id' 'transcript_id(s)' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/tpm_table.txt > ${pathout}/process/rsem_tpm_genes.txt
rm -rdf ${pathout}/process/tmp

mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,6 ${pathout}/process/${folders[$i]}/rsem.isoforms.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/tpm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/tpm_all.txt > ${pathout}/process/tmp/tpm_table.txt
headarr=('transcript_id' 'gene_id' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/tpm_table.txt > ${pathout}/process/rsem_tpm_isoforms.txt
rm -rdf ${pathout}/process/tmp

# Generate summary table
echo "=== Generating summary table... ==="
#summarize_RNAseq_RSEM.pl is added in /usr/local/bin/
summarize_RNAseq_STAR.RSEM.pl ${pathout}/process/ ${pathout}/process/summary.txt
mkdir ${pathout}/tables
mv ${pathout}/process/*txt* ${pathout}/tables/

echo "=== Finished! ==="
