# this script adds the range of barcode values to the end of the variant data file

#in_file = open("../../20-01-28_all_runs_plus_noncoding/promoter_var_analysis/1_variant_file_plus_promoter_annotations.txt","r")
in_file = open("bc_strain_allele_map_mutations_fitness_median.txt","r")
out_file = open("bc_strain_allele_map_mutations_fitness_median_range.txt", "w+")
for line in in_file:
	if "Allele" not in line:
		if "NA" not in line:
			line = line.strip()
			line2 = line.split()
			fitness_vals = line2[11].split(",")
			for i in range(len(fitness_vals)):
				fitness_vals[i] = fitness_vals[i].replace('"',"")
			fitness_vals = list(map(float,fitness_vals))
			min_val = min(fitness_vals)
			max_val = max(fitness_vals)
			range_vals = abs(max_val - min_val)
			line += "\t" + str(min_val) + "\t" + str(max_val) + "\t" + str(range_vals) + "\n"
		else:
			line = line.strip()
			line += "\tNA\tNA\tNA\n"
		out_file.write(line)
	else:
		line = line.strip()
		line += "\tMinFitVals\tMaxFitVals\tRangeFitVals\n"
		out_file.write(line)

in_file.close()
out_file.close()
		