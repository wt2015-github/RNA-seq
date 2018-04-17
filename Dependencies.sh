# Objective: Install dependencies and make them executable via command line
# Script author: David Chen (@ydavidchen)
# Notes:
# 1. This script is for MacOS. All algorithms also run on Unix/Linux; installation method/format may be different.
# 2. Depending on your system, you may need to install additional softwares/commands (e.g. Git).

#-------------------------------Phase 1: Workspace setup-------------------------------
## Set current working directory:
## Replace this line with your computer path
cd ~/repos/RNA-seq

## Create sub-directories to store files and annotations:
mkdir iGenomes
mkdir iGenomes/Homo_sapiens/
mkdir iGenomes/Homo_sapiens/UCSC/
mkdir iGenomes/Homo_sapiens/UCSC/hg19/

iGenomes/Homo_sapiens/UCSC/hg19/Annotation/
mkdir iGenomes/Homo_sapiens/UCSC/hg19/Sequence/

mkdir iGenomes/Homo_sapiens/UCSC/hg19/Sequence/RSEM_STAR_Index
mkdir iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/

## Download genomic and/or transcriptomic annotations from UCSC and Ensembl
# genome.fa
# genes.gtf
# hg19_RefSeq.bed
mv ~/Downloads/genome.fa iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/
mv ~/Downloads/genes.gtf iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
mv ~/Downloads/hg19_RefSeq.bed iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_RefSeq.bed

## Move genomic/transcriptomic annotations from Downloads folder to
cd /iGenomes/Homo_sapiens/UCSC/hg19/Sequence/RSEM_STAR_Index

#-------------------------------Phase 2: 9-step software installation, setup and testing------------------------------
## 1. fastqc:
## Download fastqc from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ into working directory
unzip ~/Downloads/fastqc_v0.11.7.zip
chmod +x FastQC/fastqc #change permission
sudo ln -s /Users/DavidKevinChen/repos/RNA-seq/FastQC/fastqc /usr/local/bin/fastqc #absolute path required!
fastqc --help #check


## 2. multqc
pip install --upgrade pip
pip install multiqc #1-step installation into command line; cython may be required
multiqc --help #check

## 3. Trim Galore!
## Install from command line
pip install cutadapt
cutadapt --version #1.16
fastqc -v
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
sudo mv TrimGalore-0.4.5/trim_galore /usr/local/bin
trim_galore --help

## 4. fastq_screen
## https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqscreen
tar xvzf ~/Downloads/fastq_screen_v0.11.4.tar.gz
chmod +x fastq_screen_v0.11.4/fastq_screen
sudo mv fastq_screen_v0.11.4/fastq_screen /usr/local/bin
fastq_screen --help

## 5. STAR
brew install cmake gcc #if already installed, run: brew upgrade gcc
cmake -DCMAKE_CXX_COMPILER=g++-6 CMakeLists.txt
make
pwd
# git clone https://github.com/alexdobin/STAR.git #or manually download if `git` not available
# cd STAR/source
# make STARforMacStatic CXX=/usr/local/bin/gcc-7
# wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
# tar -xzf 2.5.3a.tar.gz
# cd STAR-2.5.3a
# make manual #build documentation
# cd ..
conda install star -c bioconda #if above commented code fails; but requires Anaconda/miniconda
STAR --help #check

## 6. RSEM
tar xvzf ~/Downloads/RSEM-1.3.0.tar.gz
chmod +x RSEM-1.3.0/rsem-*
sudo mv RSEM-1.3.0/rsem-* /usr/local/bin
sudo mv RSEM-1.3.0/rsem_perl_utils.pm /usr/local/bin
rsem-calculate-expression --help #check

## 7. samtools
## Download source: http://www.htslib.org/download/
tar xvzf ~/Downloads/samtools-1.8.tar.bz2
cd samtools-1.8/
./configure --prefix=/Users/DavidKevinChen/samtools #creats configure file
make
make install
cd /Users/DavidKevinChen/samtools/bin
chmod +x samtools
sudo mv samtools /usr/local/bin
samtools -v #check

## 8. Picard tools from Broad Institute
## Download source code from https://broadinstitute.github.io/picard/
## github.com/broadinstitute/picard/releases/tag/2.18.2
java -version #java version "1.8.x" required
git clone https://github.com/broadinstitute/picard.git

## Build Picard:
brew install gradle
cd picard/
./gradlew shadowJar
java -jar build/libs/picard.jar

## Set environment variable (i.e. alias):
atom ~/.bashrc #other editors OK
## Add the following to the file (without hastag; and modify path!)
# export PICARD="/Users/DavidKevinChen/repos/RNA-seq/picard/build/libs/picard.jar";

##
source ~/.bashrc

## Test Picard
java -jar $PICARD -h

## 9. Make summarize_RNAseq_STAR.RSEM.pl available for execution
sudo cp summarize_RNAseq_STAR.RSEM.pl /usr/local/bin/
