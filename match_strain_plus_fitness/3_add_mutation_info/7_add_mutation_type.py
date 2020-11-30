# this script looks at mutation types

infile = open("strain_allele_map_corrected_3.txt","r")
outfile = open("strain_allele_map_corrected_4.txt","w+")
for line in infile:
	if "Allele" not in line:
		line = line.strip().split()
		# looking at coding mutations only
		if line[9] == "NA":
			line.append("Synonymous")
		elif "*" in line[9]:
			line.append("Nonsense")
		else:
			line.append("Nonsynonymous")
		# looking at promoter mutations
		if line[7] != "NA":
			#print(line[7])
			muts = line[7].replace('"','')
			#print(muts)
			muts = muts.strip().split(",")
			promoter_mutations = []
			for mut in muts:
				#print(mut)
				if "ins" in mut:
					location = int(mut[3:-1])
				elif "del" in mut:
					location = int(mut[3:-1])
				else:
					location = int(mut[1:-1])
				if location < 0:
					promoter_mutations.append(mut)
			if len(promoter_mutations) > 0:
				mutations = ",".join(promoter_mutations)
			else:
				mutations = "NA"
			line.append(mutations)
		else:
			line.append("NA")
		line = "\t".join(line) + "\n"
	else:
		line = line.strip() + "\tAllele_Type\tPromoter_Changed\n"
	outfile.write(line)

