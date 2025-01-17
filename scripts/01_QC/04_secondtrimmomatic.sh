#!/bin/bash
#SBATCH --job-name=secondtrimmomatic
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=150G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

#################################################################
# Trimmomatic
#################################################################

module load Trimmomatic/0.39
module load parallel/20180122

# set input/output directory variables
INDIR=../../data/Walnut_Colon_RNASeq
TRIMDIR=../../results/01_QC/trim_reads_2
mkdir -p $TRIMDIR

export TMPDIR=tmp
# adapters to trim out
ADAPTERS=/isg/shared/apps/Trimmomatic/0.39/adapters/NexteraPE-PE.fa

# sample list

find $INDIR/*R1* -exec basename {} _L002_R1.fastq.gz \; > samples.txt

SAMPLES=samples.txt

# run trimmomatic

cat $SAMPLES | parallel -j 10 \
java -jar $Trimmomatic PE \
        ${INDIR}/{}_L002_R1.fastq.gz \
        ${INDIR}/{}_L002_R2.fastq.gz \
        ${TRIMDIR}/{}_trim_1.fastq.gz ${TRIMDIR}/{}_trim_orphans_1.fastq.gz \
        ${TRIMDIR}/{}_trim_2.fastq.gz ${TRIMDIR}/{}_trim_orphans_2.fastq.gz \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:25 MINLEN:36
