#!/bin/bash
#SBATCH --job-name=arraytest
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 3
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[1-20]%10

hostname
date

INDIR=../../results/01_QC/trim_reads_3
OUTDIR=test
mkdir -p ${OUTDIR}

# this is an array job. 
    # one task will be spawned for each sample
    # for each task, we specify the sample as below
    # use the task ID to pull a single line, containing a single accession number from the accession list
    # then construct the file names in the call to hisat2 as below


ACCLIST=../01_QC/samples.txt

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${ACCLIST})

tail $INDIR/${SAMPLE}_trim_1.fastq.gz > $OUTDIR/${SAMPLE}

