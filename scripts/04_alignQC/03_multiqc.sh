#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Aggregate reports using MultiQC
#################################################################

module load MultiQC/1.9


qualimap=../../results/04_alignQC/qualimap
samtools=../../results/04_alignQC/samtools_stats

samtoolsmulti=../../results/04_alignQC/multiqc/samtools
mkdir -p $samtoolsmulti 
qualimapmulti=../../results/04_alignQC/multiqc/qualimap 
mkdir -p $qualimapmulti 

# run on samtool stats output
multiqc -f -o $samtoolsmulti $samtools

# run on qualimap output
multiqc -f -o $qualimapmulti $qualimap
