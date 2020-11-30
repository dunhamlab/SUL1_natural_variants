import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import re

sample = "64" # sample name 

output_61 = open("YCY"+sample+"_pool_bc_counts.txt", "w+") # format file should be in
output_61.write("barcode\tT1\tT2\tT3\tT4\tT5\tT6\tT7\tT8\tT9\tT10\tT11\n") # writes tab-delimited file with barcode and frequency of barcodes at each time point
T0_bc_counts = {}
T0_bc_counts_normalized = {}

# this counts all the samples in first time point, which is used for normalizing
if not os.path.exists("T0_"+sample+"_counts.txt"): # only run if hasn't been run before. because it takes forever!
	T0_file = open("paired_fastqs/"+sample+"_T0_merged.assembled.fastq","r")
	T0_totalReads = 0
	iter = 0
	while True:
		T0_file.readline()
		line = T0_file.readline() # this is the line with the barcode sequence. we don't really care about the rest of the lines in the file
		T0_file.readline()
		T0_file.readline()
	
		if not line:
			break
	
		line = line.split()
		line = line[0]
	
		if line in T0_bc_counts.keys():
			T0_bc_counts[line] += 1
		else:
			T0_bc_counts[line] = 1
		T0_totalReads += 1
	
		if iter % 5000 == 0:
			print(iter)
		iter += 1

	sample_output = open("T0_"+sample+"_counts.txt","w+")
	for key in T0_bc_counts:
		T0_bc_counts_normalized[key] = float(T0_bc_counts[key])/float(T0_totalReads) # frequency of barcodes
		sample_output.write(key+"\t"+str(T0_bc_counts[key])+"\t"+str(T0_bc_counts_normalized[key]) + "\t" + str(T0_totalReads)+"\n")
	T0_file.close()
else:
	T0_file = open("T0_"+sample+"_counts.txt", "r") # reads in frequency file if already done
	for line in T0_file:
		line = line.split()
		bc = line[0]
		T0_bc_counts[bc] = float(line[1])
		T0_bc_counts_normalized[bc] = float(line[2])
		T0_totalReads = float(line[3])
	print("File exists, done reading: " + str(T0_totalReads) + "\t" + str(len(T0_bc_counts_normalized.keys())))
	T0_file.close()


# counts barcodes for all other timepoints, normalizing to T0 during. creates file with all counts, each time point separate
normalized_T1_T11 = {}
for key in T0_bc_counts_normalized.keys():
	normalized_T1_T11[key] = []
rangeVal = 12
for tp in range(1,rangeVal):
	print("Time Pt: " + str(tp))
	filename = "paired_fastqs/"+sample+"_T"+str(tp)+"_merged.assembled.fastq"
	outfile_name = "T"+str(tp)+"_"+sample+"_counts.txt"
	if not os.path.exists(outfile_name):
		file = open(filename, "r")
		totalReads = 0
		bc_counts = {}
		bc_counts_normalized = {}
		while True:
			file.readline()
			line = file.readline()
			file.readline()
			file.readline()
		
			if not line:
				break
		
			line = line.split()
			line = line[0]
		
			if line in bc_counts.keys():
				bc_counts[line] += 1
			else:
				bc_counts[line] = 1
		
			if totalReads % 5000 == 0:
				print(totalReads)
			totalReads += 1
		
		sample_output = open(outfile_name,"w+")

		for key in bc_counts:
			if key in T0_bc_counts_normalized.keys():
				bc_counts_normalized[key] = (float(bc_counts[key])/float(totalReads))/T0_bc_counts_normalized[key]
				sample_output.write(key+"\t"+str(bc_counts[key])+"\t"+str(bc_counts_normalized[key]) + "\t" + str(totalReads)+"\n")
				#normalized_T1_T11[key] = []
				normalized_T1_T11[key].append(bc_counts_normalized[key])
# 			elif key in T0_bc_counts_normalized.keys() and key not in normalized_T1_T11.keys():
# 			else:
# 				normalized_T1_T11[key].append(None)
# 				sample_output.write(key+"\t"+str(bc_counts[key])+"\t"+"None"+"\t"+str(totalReads)+"\n")
		file.close()
		sample_output.close()
	else:
		print("Time Pt: " + str(tp) + " - file exists. Reading file.")
		file = open(outfile_name, "r")
		for line in file:
			line = line.split()
			bc = line[0]
# 			if line[2] == "None":
# 				line[2] = None
			normalized_T1_T11[bc].append(line[2])
			#print(normalized_T1_T11[key])
		file.close()

# merges time points into one file. 
print(len(normalized_T1_T11))
for bc in normalized_T1_T11.keys():
# 	if len(normalized_T1_T11[bc]) != rangeVal-1:
# 		print(bc + " doesn't have correct number of values.")
# 	else:
	#print(len(normalized_T1_T11[bc]))
	if len(normalized_T1_T11[bc]) == rangeVal-1:
		output_61.write(bc)
		#print(len(normalized_T1_T11[bc]))
		for i in normalized_T1_T11[bc]:
			output_61.write("\t"+str(i))
		output_61.write("\n")

output_61.close()

# for time in range(12):
# 	file = open("paired_fastqs/"+sample+"_T"+str(time)+"_merged.assembled.fastq")
	
		