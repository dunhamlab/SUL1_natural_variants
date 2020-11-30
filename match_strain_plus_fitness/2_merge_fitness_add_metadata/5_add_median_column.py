# this script takes the median instead of the average fitness of barcodes

import numpy

in_file = open("bc_strain_allele_map_plusmutations_plusfitness.txt","r")
out_file = open("bc_strain_allele_map_mutations_fitness_median.txt","w+")
for line in in_file:
	if "Allele" in line:
		line = line.strip()
		line += "\tmedian_lowfilter_fitness\n"
	else:
		line2 = line.strip().split()
		fitnesses = line2[11]
		fitnesses = fitnesses.replace('"',"")
		fitnesses = fitnesses.split(",")
		if "NA" not in fitnesses:
			fitnesses = list(map(float,fitnesses))
			#print(type(fitnesses[0]))
			median = numpy.median(fitnesses)
			line = line.strip()
			line += "\t" + str(median) + "\n"
		else:
			line = line.strip()
			line += "\tNA\n"
	out_file.write(line)

in_file.close()
out_file.close()