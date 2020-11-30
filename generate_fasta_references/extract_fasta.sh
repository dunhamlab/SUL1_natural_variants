#!/usr/bin/env bash
#$ -l mfree=4G
#$ -l h_rt=4:0:0
#$ -cwd
#$ -N get_coding_fasta
#$ -m abe
#$ -M ccyeh@uw.edu

## extract regions
module load samtools/1.9
module load tabix/latest

#use reference genome that is annotated with numbers instead of roman numerals
REF=/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-05-15_generate_CDS_with_noncoding/from_PJ/sacCer3_num.fasta
#OUTPUT_REF=SUL1_ref.fasta #change
#POS=2:788391-792075 #change. this is the region you want.
#OUTPUT=../SUL1_plusNoncoding.fasta #change

OUTPUT_REF=SUL1_CDS_ref.fasta
POS=2:789235-791814 # this needs to be changed for genomic location
OUTPUT=../SUL1_CDS.fasta

# generate reference fasta file from whole genome. output will be OUTPUT_REF
samtools faidx $REF -o ${OUTPUT_REF} ${POS}

#source loadBcftools
module load htslib/1.9
module load bcftools/1.9 

##### this extract all variants

#do tabix once only!
tabix -p from_PJ/1011Matrix_chr.vcf.gz

#touches file
> ${OUTPUT} #this is ../SUL1_plusNoncoding.fasta

#don't need to change anything here
for j in `cat strain_selected.txt`
do 
	bcftools consensus -f ${OUTPUT_REF} -I -s $j /net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-05-15_generate_CDS_with_noncoding/from_PJ/1011Matrix_chr.vcf.gz >> ${OUTPUT} #../SUL1_plusNoncoding.fasta
done

python rename_fasta.py

#rm ../SUL1_plusNoncoding.fasta

##### this only extracts SNPs
# 
# for j in `cat strain_selected.txt`
# do
# 	bcftools consensus -f SUL1_ref.fasta -s $j ../from_PJ/1011Matrix_SNPs_only_ID_chr_no_biallele.recode.vcf.gz >> SUL1_noncoding_SNPs.fasta
# done
# 
# python ../from_PJ/rename_fasta_main2.py

