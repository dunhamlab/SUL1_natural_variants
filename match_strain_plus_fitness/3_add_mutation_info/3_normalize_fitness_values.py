# this script normalizes the fitness values to wild type

import numpy

infile = open("strain_allele_map_corrected_2.txt","r")
for line in infile:
	line = line.strip().split()
	if line[1] == "863":
		average_WT_fitness = float(line[12])
infile.close()
#average_WT_fitness = 0.009699439

infile = open("strain_allele_map_corrected_2.txt","r")
outfile = open("strain_allele_map_corrected_3.txt","w+")
for line in infile:
	#line = line.replace('"','')
	line = line.strip().split()
	stdev = "NA"
	if line[0] != "Allele" and line[11] != "NA":
		bc_list = line[11]
		bc_list = bc_list.replace('"','')
		bc_list = bc_list.split(",")
		print(bc_list)
		for i in range(len(bc_list)):
			bc_list[i] = float(bc_list[i])
			bc_list[i] = bc_list[i] - average_WT_fitness
		print(bc_list)
		average = numpy.mean(bc_list)
		median = numpy.median(bc_list)
		stdev = numpy.std(bc_list)
		minfit = min(bc_list)
		maxfit = max(bc_list)
		rangefit = maxfit-minfit
		bc_list = map(str,bc_list)
		line[11] = ",".join(bc_list)
		line[12] = str(average)
		line[13] = str(median)
		line[14] = str(minfit)
		line[15] = str(maxfit)
		line[16] = str(rangefit)
	elif line[0] == "Allele":
		stdev = "StDev"
	line = "\t".join(line) + "\t" + str(stdev) + "\n"
	outfile.write(line)

outfile.close()