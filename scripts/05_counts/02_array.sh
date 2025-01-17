#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[0-15]%16

hostname
date

#################################################################
# Generate Counts 
#################################################################
module load htseq/0.13.5


INDIR=../../results/03_align/alignments
OUTDIR=../../results/05_counts/counts
mkdir -p $OUTDIR

# accession list

SAMPLES=../01_QC/samples.txt
NUM=$(expr ${SLURM_ARRAY_TASK_ID} + 1)

SAMPLE=$(sed -n ${NUM}p $SAMPLES)

# gtf formatted annotation file
GTF=../../data/genome/Homo_sapiens.GRCh38.112.gtf

# run htseq-count on each sample in parallel
htseq-count \
        -s reverse \
        -r pos \
        -f bam $INDIR/${SAMPLE}.bam \
        $GTF \
        > $OUTDIR/${SAMPLE}.counts
