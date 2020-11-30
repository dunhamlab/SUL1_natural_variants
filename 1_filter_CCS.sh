#$ -cwd
#$ -N CCS_analysis
#$ -o output_map.txt
#$ -e error_map.txt
#$ -m abe
#$ -M ccyeh@uw.edu

module load samtools/latest

#create fasta file, filtering for a certain number of passes. example shown below is for minimum 4 passes. the number in {print $13} can be 12 or 13. so check by:
#this is the CCS bam file automatically (hopefully) run by the Eichler lab. I don't have the scripts for how they generated these.
#CHANGE
# BAM=/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/data/181002_PacBio_Sequencing_results/YehC/SUL1_bc_allele_lib.CY/CCS_001279/SUL1_bc_allele_lib.CY.ccs.bam
BAM=/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/data/190820_PacBio_Sequencing_results_2/SUL1_bc_allele_lib_2_CY-6.CY/CCS_002444/m54179_190810_110541.ccs.bam

samtools view ${BAM} | awk '{print $13}' | less
#once you've checked that, create a new fastq file that filters for reads with at least 10 passes. The number 5 in the command below indicates the minimum number of passes. Change that, and the output file, as appropriate.
samtools view ${BAM} | awk '{split($13,a,":");if (a[3]>2) print}' | samtools fastq - > ccs_np3.fastq
