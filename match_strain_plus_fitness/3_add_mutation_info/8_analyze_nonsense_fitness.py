# this script analyzes nonsense mutations and whether or not location of it affects function.

import numpy

infile = open("strain_allele_map_corrected_4.txt","r")

stop_fitness_dict = {}
for line in infile:
	line = line.strip().split("\t")
	if line[11] != "bc_fitnesses" and line[11] != "NA":
		if "*" in line[9]:
			mutations = line[9]
			mutations = mutations.split(",")
			for mut in mutations:
				if "*" in mut:
					mut = mut.replace('"','')
					print(line[0] + "\t" + mut)
					loc = int(mut[1:-1])
					print(loc)
					break
			fitnesses = line[11]
			if loc not in stop_fitness_dict:
				stop_fitness_dict[loc] = fitnesses
			else:
				stop_fitness_dict[loc] = stop_fitness_dict[loc] + "," + fitnesses

stop_fitness_mean = {}
stop_fitness_median = {}
stop_fitness_stdev = {}
for loc in stop_fitness_dict:
	#print(loc)
	print(str(loc) + "\t" + stop_fitness_dict[loc])
	fitness_vals = stop_fitness_dict[loc].split(",")
	for i in range(len(fitness_vals)):
			fitness_vals[i] = float(fitness_vals[i])
	stop_fitness_mean[loc] = numpy.mean(fitness_vals)
	stop_fitness_median[loc] = numpy.median(fitness_vals)
	stop_fitness_stdev[loc] = numpy.std(fitness_vals)
	
outfile = open("stop_location_fitness_fixed.txt","w+")
for loc in stop_fitness_mean:
		#print(str(loc) + "\t" + str(stop_fitness_mean[loc]) + "\n")
		outfile.write(str(loc) + "\t" + str(stop_fitness_dict[loc]) + "\t" + str(stop_fitness_mean[loc]) + "\t" + str(stop_fitness_median[loc]) + "\t" + str(stop_fitness_stdev[loc]) + "\n")
