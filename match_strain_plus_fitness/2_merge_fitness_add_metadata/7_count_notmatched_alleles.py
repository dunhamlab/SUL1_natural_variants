# this script counts the number of unmatched alleles in pacbio, and how many barcodes are associated with those

infile = open("strain_bc_map.tsv","r")
notmatched = 0
matched = 0
for line in infile:
	line = line.strip().split("\t")
	#print(line)
	if "Allele" not in line:
		#print(line[5])
		if line[5] == "":
			notmatched += int(line[1])
			print(line[1])
		else:
			matched += int(line[1])
			
print(notmatched)
print(matched)

not_matched_frac = notmatched/(notmatched+matched)
print(not_matched_frac)
			
			
		
