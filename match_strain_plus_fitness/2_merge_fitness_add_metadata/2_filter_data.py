# this file removes barcodes that aren't in the three replicates

dataset = open("strain_bc_map.tsv","r")
fitness_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-08-13_barseq_60_62/9_191126_merged_slopes_0.3_cutoff_no60_filtered.tsv","r")

new_fitness = open("updated_fitness.txt","w+")
valid_bcs = []
new_fitness.write("barcode\t61_slope\t64_slope\t62_slope\taverage\n")
fitness_line_dict = {}
for line in fitness_file:
	line = line.strip().split()
	if line[0] != "barcode":
		#if float(line[5]) < 0.08:
		valid_bcs.append(line[0])
		new_fitness.write("\t".join(line)+"\n")
		fitness_line_dict[line[0]] = "\t".join(line)+"\n"
fitness_file.close()
new_fitness.close()

new_dataset = open("updated_strain_bc_map.tsv","w+")
#new_dataset.write("Allele\tNum_BCs\tAllele_Len\tAllele\tBCs\tStrain\tNoAmbiguity")
all_lines = []
all_bcs = []
for line in dataset:
	line = line.strip().split()
	if line[0] != "Allele":
		bcs_w_fitness = []
		bcs = line[4].split(",")
		for bc in bcs:
			if bc in valid_bcs:
				bcs_w_fitness.append(bc)
		line[4] = ",".join(bcs_w_fitness)
	if len(line) == 7:
		new_dataset.write("\n"+"\t".join(line))
		new_bcs = line[4].split(",")
		all_bcs = all_bcs + new_bcs
# 	else:
# 		new_dataset.write("\n"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t\t"+line[5])
new_dataset.close()

#filter dataset for ones that are also in haplotype map
all_bcs = list(set(all_bcs))
for item in all_bcs:
	print(item)
print(len(all_bcs))
# fitness_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-08-13_barseq_60_62/9_191126_merged_slopes_0.3_cutoff_no60_filtered.tsv","r")
#new_dataset = open("191126_updated_strain_bc_map.tsv","r")
new_fitness_2 = open("fitness_in_pb_results.txt","w+")
new_fitness_2.write("barcode\t61_slope\t64_slope\t62_slope\taverage\tstdev\n")
new_fitness = open("updated_fitness.txt","r")
for line in new_fitness:
	line2 = line.strip().split()
	if line2[0] in all_bcs:
		new_fitness_2.write(line)
new_fitness.close()
new_fitness_2.close()

# removed some psc alleles
# new_dataset = open("updated_strain_bc_map.tsv","r")
# pscs_bcs = []
# for line in new_dataset:
# 	line2 = line.strip().split()
# 	if len(line2) > 0:
# 		#print(line2)
# 		if line2[0] == "190" or line2[0] == "95" or line2[0] == "260" or line2[0] == "369":
# 			bcs = line2[4].split(",")
# 			pscs_bcs = pscs_bcs + bcs
# 
# no_psc_fitness = open("fitness_nopsc_bcs.txt","w+")
# new_fitness_2 = open("fitness_in_pb_results.txt","r")
# for line in new_fitness_2:
# 	line2 = line.strip().split()
# 	if line2[0] not in pscs_bcs:
# 		no_psc_fitness.write(line)
# no_psc_fitness.close()
new_fitness_2.close()
	