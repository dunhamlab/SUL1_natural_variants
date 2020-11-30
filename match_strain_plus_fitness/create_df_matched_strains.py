# this script creates a dataframe of strains that matched with the correct/expected number of alleles. (haploid = 1 match, diploid = <3 matches, etc.)
# metadata should be added here too?

from metadata_analysis import strain_metadata

haploid_infile = open("output/haploid_matches.txt","r")
haploid_dict = {}
allele_dict_hap = {}
for line in haploid_infile:
	line = line.strip().split()
	strains_str = line[0]
	strains = line[0].split(",")
	seq = line[2]
	for strn in strains:
		haploid_dict[strn] = seq
	allele_dict_hap[seq] = strains

# print("haploids")
# for item in allele_dict_hap:
# 	print(allele_dict_hap[item])
	

diploid_infile = open("output/diploid_matches.txt","r")
diploid_dict = {}
allele_dict_dip = {}

lineread = diploid_infile.readline()
while True:
	strain = lineread[1:].strip()
	seqs = []
	lineread = diploid_infile.readline()
	while ">" not in lineread:
		seqs.append(lineread.strip())
		lineread = diploid_infile.readline()
		if lineread == "":
			break
	if lineread == "":
		break
	if len(seqs) < 3:
		diploid_dict[strain] = seqs
		for allele in diploid_dict[strain]:
			if allele not in allele_dict_dip:
				allele_dict_dip[allele] = [strain]
			else:
				allele_dict_dip[allele].append(strain)

# print("diploids")
# for item in allele_dict_dip:
# 	print(allele_dict_dip[item])
# 	
			
triploid_infile = open("output/triploid_matches.txt","r")
triploid_dict = {}
allele_dict_trip = {}

lineread = triploid_infile.readline()
while True:
	strain = lineread[1:].strip()
	seqs = []
	lineread = triploid_infile.readline()
	while ">" not in lineread:
		seqs.append(lineread.strip())
		lineread = triploid_infile.readline()
		if lineread == "":
			break
	if lineread == "":
		break
	if len(seqs) < 4:
		triploid_dict[strain] = seqs
		for allele in triploid_dict[strain]:
			if allele not in allele_dict_trip:
				allele_dict_trip[allele] = [strain]
			else:
				allele_dict_trip[allele].append(strain)			

# print("triploids")
# for item in allele_dict_trip:
# 	print(allele_dict_trip[item])

tetraploid_infile = open("output/tetraploid_matches.txt","r")
tetraploid_dict = {}
allele_dict_tetra = {}

lineread = tetraploid_infile.readline()
while True:
	strain = lineread[1:].strip()
	seqs = []
	lineread = tetraploid_infile.readline()
	while ">" not in lineread:
		seqs.append(lineread.strip())
		lineread = tetraploid_infile.readline()
		if lineread == "":
			break
	if lineread == "":
		break
	if len(seqs) < 5:
		tetraploid_dict[strain] = seqs
		for allele in tetraploid_dict[strain]:
			if allele not in allele_dict_tetra:
				allele_dict_tetra[allele] = [strain]
			else:
				allele_dict_tetra[allele].append(strain)		

# print("tetraploids")
# for item in allele_dict_tetra:
# 	print(allele_dict_tetra[item])			
							
pentaploid_infile = open("output/pentaploid_matches.txt","r")
pentaploid_dict = {}
allele_dict_pent = {}

lineread = pentaploid_infile.readline()
while True:
	strain = lineread[1:].strip()
	seqs = []
	lineread = pentaploid_infile.readline()
	while ">" not in lineread:
		seqs.append(lineread.strip())
		lineread = pentaploid_infile.readline()
		if lineread == "":
			break
	if lineread == "":
		break
	if len(seqs) < 5:
		pentaploid_dict[strain] = seqs
		for allele in pentaploid_dict[strain]:
			if allele not in allele_dict_pent:
				allele_dict_pent[allele] = [strain]
			else:
				allele_dict_pent[allele].append(strain)		

# print("pentaploids")
# for item in allele_dict_pent:
# 	print(allele_dict_pent[item])		

print("unique alleles in haploids: " + str(len(allele_dict_hap)))
print("unique alleles in diploids: " + str(len(allele_dict_dip)))
print("unique alleles in triploids: " + str(len(allele_dict_trip)))
print("unique alleles in tetraploids: " + str(len(allele_dict_tetra)))
print("unique alleles in pentaploids: " + str(len(allele_dict_pent)))

all_alleles = list(allele_dict_hap.keys()) + list(allele_dict_dip.keys()) + list(allele_dict_trip.keys()) + list(allele_dict_tetra.keys()) + list(allele_dict_pent.keys())
set_alleles = list(set(all_alleles))
print("unique alleles combined: " + str(len(set_alleles)))

list1 = [allele_dict_hap,allele_dict_dip,allele_dict_trip,allele_dict_tetra,allele_dict_pent]

# for dict in list1:
# 	for allele in dict:
# 		print(dict[allele])

# IDK why i decided to recreate this. i already have it!! ugh
allele_df = open("allele_df.txt","w+")
allele_df.write("allele_identifier\tstrains\tallele")
num = 1
for allele in set_alleles:
	strains = []
	print(num)
	if allele in allele_dict_hap:
		#print(allele_dict_hap[allele][0])
		strains.append(allele_dict_hap[allele][0])
	if allele in allele_dict_dip:
		if len(allele_dict_dip[allele]) > 1:
			for item in allele_dict_dip[allele]:
				strains.append(item)
		else:
			strains.append(allele_dict_dip[allele][0])
	if allele in allele_dict_trip:
		if len(allele_dict_trip[allele]) > 1:
			for item in allele_dict_trip[allele]:
				strains.append(item)
		else:
			strains.append(allele_dict_trip[allele][0])	
	if allele in allele_dict_tetra:
		if len(allele_dict_tetra[allele]) > 1:
			for item in allele_dict_tetra[allele]:
				strains.append(item)
		else:
			strains.append(allele_dict_tetra[allele][0])
	if allele in allele_dict_pent:
		if len(allele_dict_pent[allele]) > 1:
			for item in allele_dict_pent[allele]:
				strains.append(item)
		else:
			strains.append(allele_dict_pent[allele][0])
	strains = ",".join(strains)
	allele_df.write("\n"+str(num)+"\t"+strains+"\t"+allele)
	num += 1
# allele_df.close()







