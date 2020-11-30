# this script looks at variance across alleles

import numpy

infile = open("strain_allele_map_corrected_4.txt","r")

allele_variance_dict = {}
allele_mean_dict = {}
for line in infile:
	line = line.strip().split()
	if line[0] != "Allele":
		line[16] = line[16].replace('"','')
		print(line[16])
		print(type(line[16]))
		if line[16] != "0.0":
			if line[16] != "NA":
				fitnesses = line[11].split(",")
				for i in range(len(fitnesses)):
					fitnesses[i] = float(fitnesses[i])
				variance = numpy.std(fitnesses)
				allele_mean_dict[line[0]] = float(line[12])
				allele_variance_dict[line[0]] = variance

outfile = open("allele_bc_variance.txt","w+")
for key in allele_mean_dict:
	if allele_mean_dict[key] <= -0.07:
		outfile.write(key+"\t"+str(allele_variance_dict[key])+"\tlof\n")
	elif allele_mean_dict[key] <= -0.02 and allele_mean_dict[key] > -0.07:
		outfile.write(key+"\t"+str(allele_variance_dict[key])+"\tpromoterchange\n")
	else:
		outfile.write(key+"\t"+str(allele_variance_dict[key])+"\twt\n")
		