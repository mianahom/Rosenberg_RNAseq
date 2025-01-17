#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=100G
#SBATCH --partition=general
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

module load samtools/1.16.1
module load parallel/20180122

export TMPDIR=tmp

ALIGNMENTS=../../results/03_align/alignments


SAMPLES=../01_QC/samples.txt

cat $SAMPLES | parallel -j 6 \
"samtools view -@ 1 -S -h -u $ALIGNMENTS/{}.sam | samtools sort -@ 1 -T {} > $ALIGNMENTS/{}.bam"

