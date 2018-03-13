#!/bin/bash
#Ting Wang, 201606

manual="
Shell script for calling variants from results of RNA-seq analysis, after running run_RNAseq.sh
Usage:
    callSNP_RNAseq.sh -g <hg19> -i <input_folder> -t <threads>
Options:
    -g    default is hg19, only support hg19 now because vcf files of other genomes are not available
    -i    directory of inputs, its the resultant folder of run_RNAseq.sh
    -t    average number of threads for each sample, must be integer, default is 1
"

if [[ $# -le 0 ]]; then
    echo "${manual}" >&2
    exit 1
fi

while getopts :g:i:t: ARGS
do
case $ARGS in
    g)
        genome="hg19"
        ;;
    i)
        pathin=$OPTARG
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

echo ${genome:="hg19"} >/dev/null
echo ${pathin:="nopathin"} >/dev/null
echo ${threads:=1} >/dev/null

if [[ ${pathin} == "nopathin" ]]; then
    echo "=== please set directory of inputs!! === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathin} ]]; then
    echo "=== input directory does not exist!! === ${manual}" >&2
    exit 1
fi

(( ${threads} )) 2>/dev/null
if [[ $? != 0 || ! ${threads} -ge 1 ]]; then
    echo "=== thread number should be a positive integer!! === ${manual}" >&2
    exit 1
fi
echo ""
echo "=== set ${threads} threads for each sample ==="

if [[ ${genome} == "hg19" ]]; then
    echo "=== set genome: hg19, other genomes are not supported now, because vcf files are unavailable ==="
    genome_fa="/home1/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    dbsnp_vcf="/home1/data/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/GATK_hg19bundle/dbsnp_138.hg19.vcf"
    indel_vcf="/home1/data/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/GATK_hg19bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
else
    echo "=== only support hg19 now, because vcf files of other genomes are not available === ${manual}" >&2
    exit 1
fi

folders=(`ls -d ${pathin}/*/ | sort`)
if [[ ${#folders[@]} -eq 0 ]]; then
    echo "=== there is no subfolder in the directory!! === ${manual}" >&2
    exit 1
else
    echo "=== there are ${#folders[@]} RNAseq resultant subfolders in the directory ==="
fi

for ((i=0; i<${#folders[@]}; i++)); do
    subfolder=$(echo ${folders[$i]} | rev | cut -d '/' -f2 | rev)
    echo "=== calling vairants for sample ${subfolder} ==="
    {
        mkdir ${pathin}/${subfolder}/tmp_java
        #echo "=== marking duplicates for sample ${subfolder} ==="
        #java -Xmx10g -XX:ParallelGCThreads=${threads} -jar $PICARD MarkDuplicates I=${pathin}/${subfolder}/Aligned.sortedByCoord.out.bam O=${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam M=${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=${pathin}/${subfolder}/tmp_java/
	echo "=== adding read group information to bam ==="
        java -jar $PICARD AddOrReplaceReadGroups I=${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam O=${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.bam RGID=${subfolder} RGLB=${subfolder} RGPL=illumina RGPU=${subfolder} RGSM=${subfolder} TMP_DIR=${pathin}/${subfolder}/tmp_java/
        echo "=== generate index file ==="
        samtools index ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.bam
        echo "=== splitting 'N' trim ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T SplitNCigarReads -R ${genome_fa} -I ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.bam -o ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
        echo "=== running base quality score recalibration (BQSR) 1: building recalibration model ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T BaseRecalibrator -R ${genome_fa} -I ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.bam -knownSites ${dbsnp_vcf} -knownSites ${indel_vcf} -o ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.table
        echo "=== running BQSR 2: analyzing covariation remaining after recalibration ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T BaseRecalibrator -R ${genome_fa} -I ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.bam -knownSites ${dbsnp_vcf} -knownSites ${indel_vcf} -BQSR ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.table -o ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.table.post
        echo "=== running BQSR 3: generating recalibration plots ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T AnalyzeCovariates -R ${genome_fa} -before ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.table -after ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.table.post -plots ${pathin}/${subfolder}/recalibration_plots.pdf
        echo "=== running BQSR 4: applying recalibration to seq data ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T PrintReads -R ${genome_fa} -I ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.bam -BQSR ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.table -o ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.reads.bam
        echo "=== calling variants ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T HaplotypeCaller -R ${genome_fa} -I ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.reads.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.reads.output.vcf
        echo "=== filtering variants ==="
        java -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${pathin}/${subfolder}/tmp_java/ -jar $GATK -T VariantFiltration -R ${genome_fa} -V ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.reads.output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${pathin}/${subfolder}/Aligned.sortedByCoord.out.dedupped.addRG.split.recal.reads.output.filted.vcf
        rm -rdf ${pathin}/${subfolder}/tmp_java
    } > ${pathin}/${subfolder}/log.callSNP.txt 2>&1 &
done
wait
echo "=== SNP calling Finished ==="

