#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=16G -l h_rt=1:0:0:0
#!/usr/bin/env bash
#$ -cwd
#$ -N CCS_analysis
#$ -o output_map.txt
#$ -e error_map.txt
#$ -m abe
#$ -M ccyeh@uw.edu

export PACBIO=/net/shendure/vol10/projects/subassemblyByPacBio/src/smrtanalysis/smrtcmds/bin/

module load python/2.7.3
module load pysam/0.7.5
module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1
module load bwa/0.7.13

#CHANGE. this one is the fastq that's been filtered for minimum # of passes
FASTQ=$1 #full path to input ccs fastq file (can be gzipped)
#CHANGE
REF=$2 #Full path to reference fasta
#CHANGE
OUTDIR=$3 #path to output directory
#CHANGE
SAMPLE=$4 #name of sample eg: CYP2C9_PB1

#COMMENT OR UNCOMMENT!!!
# Index reference (only run once!!)
(>&2 echo ***CreateIndexFiles***)
bwa index -a is $REF
$PACBIO/samtools faidx $REF

# Align reads:
-C (append fasta comment to output ?)  -M (mark shorter split hits as secondary) -L 80 (increase clipping penalty)
(>&2 echo ***BWA - Mem***)
/net/shendure/vol1/home/mkircher/bin/bwa_new mem -C -M -L 80 $REF $FASTQ | /net/shendure/vol1/home/mkircher/bin/samtools view -uS - | /net/shendure/vol1/home/mkircher/bin/samtools sort - ${OUTDIR}/${SAMPLE}_ccs_aligned

# # Not sure why using modified samtools for flagstat, but this is what was in script
(>&2 echo ***Samtools - Flagstat***)
/net/shendure/vol1/home/mkircher/bin/mod_samtools/samtools flagstat ${OUTDIR}/${SAMPLE}_ccs_aligned.bam > ${OUTDIR}/${SAMPLE}_ccs_aligned_stats.txt

# # Get rid of extra flags in bam file
(>&2 echo ***Samtools - View***)
$PACBIO/samtools view -h ${OUTDIR}/${SAMPLE}_ccs_aligned.bam | cut -f -14 | $PACBIO/samtools view -Sb - > ${OUTDIR}/${SAMPLE}_ccs_aligned.fix.bam

# # Check CIGAR strings (samtools version not specified here, using PacBio version)
(>&2 echo ***Samtools - CigarString***)
$PACBIO/samtools view ${OUTDIR}/${SAMPLE}_ccs_aligned.fix.bam | cut -f 6 | sort| uniq -c | sort -nr > ${OUTDIR}/${SAMPLE}_cigars.txt

#Make file of BC/sequence
#-s flag considers soft clipped reads (default is off)
#Using Melissa's edited version, but commented out the only print statement so the output can go right to a text file. 
# NEED TO INPUT COORDINATES of bc length, position, and CDS start and end **********BE SURE TO CHANGE FOR LIBRARY HERE****************
(>&2 echo ***BarcodeInsertPairs***)
# with noncoding:
/net/dunham/vol2/Clara/projects/CYP_project/PacBio/scripts/extractBarcodeInsertPairs_moreQC.py --verbose ${OUTDIR}/${SAMPLE}_ccs_aligned.fix.bam -l 10 -p 109 --start 110 --end 3795 -b 2 > ${OUTDIR}/${SAMPLE}_seq_barcodes_filtered.txt

#CDS only
#/net/dunham/vol2/Clara/projects/CYP_project/PacBio/scripts/extractBarcodeInsertPairs_moreQC.py --verbose ${OUTDIR}/${SAMPLE}_ccs_aligned.fix.bam -l 10 -p 109 --start 371 --end 2951 -b 2 > ${OUTDIR}/${SAMPLE}_seq_barcodes_filtered.txt

#unify barcodes: take most common sequence or best quality
(>&2 echo ***UnifyAssignment***)
/net/shendure/vol1/home/mkircher/bin/PacBioBarcodeAssignments/unifyAssignment.py ${OUTDIR}/${SAMPLE}_seq_barcodes_filtered.txt | gzip -c > ${OUTDIR}/${SAMPLE}_combined_minQ0_assignment.tsv.gz
gunzip ${OUTDIR}/${SAMPLE}_combined_minQ0_assignment.tsv.gz

# #From there I use that minQ0_assignment to generate a fake fastq file that I can then input into enrich2. The minQ0_assignment file itself serves as your barcode-variant map that you use in enrich2.
(>&2 echo ***MakeFakeFastQ***)
python /net/dunham/vol2/Clara/projects/CYP_project/PacBio/scripts/make_fake_fastq.py ${OUTDIR}/${SAMPLE}_combined_minQ0_assignment.tsv ${OUTDIR}/${SAMPLE}_BC_only.fastq ${SAMPLE}
python /net/dunham/vol2/Clara/projects/CYP_project/PacBio/scripts/make_fake_fastq.py ${OUTDIR}/${SAMPLE}_seq_barcodes_filtered.txt ${OUTDIR}/${SAMPLE}_BC_only_uncompressed.fastq ${SAMPLE}

# (>&2 echo ***Compress***)
gzip ${OUTDIR}/${SAMPLE}_BC_only.fastq
gzip ${OUTDIR}/${SAMPLE}_BC_only_uncompressed.fastq
