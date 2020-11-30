# this script analyzes the results by strain
# how many strains are matched?
# how many alleles per strain?
# how many barcodes per allele?

infile = open("strain_allele_map_corrected_3.txt","r")
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
outfile = open("alleles_per_strain.txt","w+")
for strain in strain_alleles:
	print(strain + "\t" + ",".join(strain_alleles[strain]))
	outfile.write(strain+"\t"+",".join(strain_alleles[strain])+"\n")
print(str(len(all_strains)) + " strains.")
print(str(bc) + " barcodes.")

