#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=30G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Index the Genome
#################################################################

# load software
module load hisat2/2.2.1

# input/output directories
OUTDIR=../../results/02_index/hisat2_index
mkdir -p $OUTDIR

GENOME=../../data/genome/Homo_sapiens.GRCh38.dna.toplevel.fa

hisat2-build -p 16 $GENOME $OUTDIR/Hsapiens

