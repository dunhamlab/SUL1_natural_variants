#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=12G -l h_rt=48:0:0
#$ -cwd
#$ -N merge_reads

#This script trims and merges paired-end reads 

module load pear/latest

cd /net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-30_barseq_fitness

FREAD=$1 # read 1, or forward read
RREAD=$2 # read 2, or reverse read
BCSIZE=$3 # size of barcode
SAMPLENAME=$4

python trim_reads.py ${FREAD} ${BCSIZE} # trim forward reads
python trim_reads.py ${RREAD} ${BCSIZE} #output is ${RREAD}_trimmed.fastq trims reverse reads

# generalize formatting
FREAD2=$(echo /net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-30_barseq_fitness/merged_fastqs/${FREAD} | cut -f 1 -d '.')
RREAD2=$(echo /net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-30_barseq_fitness/merged_fastqs/${RREAD} | cut -f 1 -d '.')

pear -v ${BCSIZE} -n ${BCSIZE} -m ${BCSIZE} -f ${FREAD2}_trimmed.fastq -r ${RREAD2}_trimmed.fastq -o ${SAMPLENAME}_merged # merges using quality scores

# get rid of extraneous files
rm ${SAMPLENAME}_merged.discarded.fastq
rm ${SAMPLENAME}_merged.unassembled.forward.fastq
rm ${SAMPLENAME}_merged.unassembled.reverse.fastq
