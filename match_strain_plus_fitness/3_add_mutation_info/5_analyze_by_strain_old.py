# this script compares the subassembly from 1/28 to this folder's

#oldfile = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/20-05-22_a_redo_allelesearch/3_allele_search_2/output/200522_combined_np3_strain_muts_map_pluslessstringentfitness.txt","r")

infile = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/20-05-22_a_redo_allelesearch/3_allele_search_2/output/200522_combined_np3_strain_muts_map_pluslessstringentfitness.txt","r")
all_strains = []
bc = 0 
strain_alleles = {} #key in strain, value is allele identifier (first column)
for line in infile:
	line = line.strip().split()
	if "Allele" not in line:
		strains = line[5]
		strains = strains.replace('"','')
		strains = strains.split(",")
		for strain in strains:
			if strain not in strain_alleles:
				strain_alleles[strain] = [line[0]]
			else:
				strain_alleles[strain].append(line[0])
		all_strains = all_strains + strains
		num = int(line[1])
		bc += num

#print(len(all_strains))
all_strains = list(set(all_strains))
outfile = open("alleles_per_strain_oldrun.txt","w+")
for strain in strain_alleles:
	print(strain + "\t" + ",".join(strain_alleles[strain]))
	outfile.write(strain+"\t"+",".join(strain_alleles[strain])+"\n")
print(str(len(all_strains)) + " strains.")
print(str(bc) + " barcodes.")





# 
# all_strains = []
# for line in oldfile:
# 	line = line.strip().split()
# 	if "Allele" not in line:
# 		strains = line[5]
# 		strains = strains.replace('"','')
# 		strains = strains.split(",")
# 		all_strains = all_strains + strains
# 
# all_strains = list(set(all_strains))
# print(str(len(all_strains)))