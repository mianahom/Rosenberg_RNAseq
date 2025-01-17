#!/bin/bash
#SBATCH --job-name=mem_intensive_htseq_count
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Generate Counts 
#################################################################
module load htseq/0.13.5
module load parallel/20180122

INDIR=../../results/03_align/alignments
OUTDIR=../../results/05_counts/counts
mkdir -p $OUTDIR

# accession list

SAMPLES=../01_QC/samples.txt

# gtf formatted annotation file
GTF=../../data/genome/Homo_sapiens.GRCh38.112.gtf

# run htseq-count on each sample in parallel
cat $SAMPLES | parallel -j 8 \
    "htseq-count \
        -s reverse \
        -r pos \
        -f bam $INDIR/{}.bam \
        $GTF \
        > $OUTDIR/{}.counts"
