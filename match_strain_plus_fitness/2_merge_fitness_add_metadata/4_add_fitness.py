# this file adds average fitness values of barcodes (and the average of each allele fitness)

import numpy 

#data_file = open("200109_bc_strain_allele_map_plusmutations.txt","r")
data_file=open("bc_strain_allele_map_plusmutations.txt","r")
fitness_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-08-13_barseq_60_62/9_191126_merged_slopes_0.3_cutoff_no60_filtered.tsv","r")
bc_fit_dict = {}
for line in fitness_file:
	line = line.strip().split()
	if line[0] != "barcode":
		bc_fit_dict[line[0]] = line[4]
fitness_file.close()
	

#out_file = open("200109_bc_strain_allele_muts_map_plusfitness.txt","w+")
out_file = open("bc_strain_allele_map_plusmutations_plusfitness.txt","w+")
for line in data_file:
	line = line.strip().split()
	if line[0] != "Allele":
		bcs = line[4]
		if bcs != "NA":
			fitnesses = []
			fitnesses_str = []
			bcs = bcs.split(",")
			for bc in bcs:
				fitnesses.append(float(bc_fit_dict[bc]))
				fitnesses_str.append(bc_fit_dict[bc])
			average_fit = numpy.mean(fitnesses)
			out_file.write("\t".join(line)+"\t"+",".join(fitnesses_str)+"\t"+str(average_fit)+"\n")
		else:
			out_file.write("\t".join(line)+"\t"+"NA"+"\t"+"NA"+"\n")
	else:
		out_file.write("\t".join(line)+"\tbc_fitnesses\taverage_fitness\n")
				
				
		