#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



echo `hostname`

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastqc/0.12.1
module load parallel/20180122
module load MultiQC/1.15

INDIR=../../data/Walnut_Colon_RNASeq
REPORTDIR=../../results/01_QC/raw/fastqc
mkdir -p $REPORTDIR
MULTIQC=../../results/01_QC/raw/multiqc
mkdir -p $MULTIQC


# run fastp in parallel, 6 samples at a time
ls $INDIR | parallel -j 6 \
    fastqc --outdir $REPORTDIR $INDIR/{}



# run on fastqc output
multiqc -f -o $MULTIQC $REPORTDIR
