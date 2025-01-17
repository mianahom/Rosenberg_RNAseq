#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=120G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.16.1
module load parallel/20180122

export TMPDIR=tmp

INDIR=../../results/01_QC/trim_reads_3
OUTDIR=../../results/03_align/alignments
mkdir -p ${OUTDIR}

INDEX=../../results/02_index/hisat2_index/Hsapiens

SAMPLES=../01_QC/samples.txt

cat $SAMPLES | parallel -j 6 \
hisat2 \
    -x ${INDEX} \
    -1 ${INDIR}/{}_trim_1.fastq.gz \
    -2 ${INDIR}/{}_trim_2.fastq.gz \
    -S ${OUTDIR}/{}.sam
