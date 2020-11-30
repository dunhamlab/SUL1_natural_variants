#this script shortens fastq files while retaining quality scores/read identifiers

import sys
import re
import random

fastq_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-30_barseq_fitness/merged_fastqs/" + sys.argv[1],"r") # replace path with directory fastq files are in
read_length = int(sys.argv[2]) # length of barcode. must be at beginning of read

basename = re.split("\.",sys.argv[1]) # get first part of filename as basename for naming files
basename = basename[0]
fastq_trimmed_name = basename+"_trimmed.fastq" # shortens read
fastq_trimmed = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-30_barseq_fitness/merged_fastqs/" + fastq_trimmed_name,"w+")

iter = 0
# bases = ["A","C","G","T"] # not needed anymore. was used for QC

while True:		# reads fastq file
	readname = fastq_file.readline() # name of read
	sequence = fastq_file.readline() # sequence to trim
	fastq_file.readline() 
	quality = fastq_file.readline()
	
	#preseq = random.choice(bases)+random.choice(bases)+random.choice(bases)+random.choice(bases)+random.choice(bases)
	if not sequence:
		break

	# keep sequence and quality scores
	sequence = sequence[0:read_length]
	quality = quality[0:read_length]
	
	# make into new fastq files
	fastq_trimmed.write(readname+sequence+"\n+\n"+quality+"\n")
	
	if iter % 100000 == 0:
		print(iter)
	iter += 1
	
fastq_file.close()
fastq_trimmed.close()

