# this script merges haplotypes after running label_haplotypes.py and manually fixing things (haplotypes_to_realign_fixed_manually.txt)

haplotypes_file = open("haplotypes_filtered.txt","r")

# for each line:
# [allele label, number of bcs, allele length, allele, barcodes]
haplotype_lines = {}
for line in haplotypes_file:
	line = line.strip().split()
	if line[0] != "Num":
		haplotype_lines[line[0]] = line[1:]

num_haps = len(haplotype_lines)
print("Num haplotypes: " + str(num_haps))
# WT is #369
num_haps += 1

#fixed_haplotypes_file = open("haplotypes_to_realign_fixed_manually.txt","r")
fixed_haplotypes_file = haplotypes_file

#fixed_haplotypes_file = open("haplotypes_to_realign.txt","r")
for line in fixed_haplotypes_file:
	line = line.strip().split()
	if line[0] != "Real_allele":
		if line[0] in haplotype_lines:
			print(len(haplotype_lines[line[0]]))
			haplotype_lines[line[0]][3] += ","+line[5]
			bc_list = haplotype_lines[line[0]][3].split(",")
			haplotype_lines[line[0]][0] = len(bc_list)
			#haplotype_lines[line[0]][0] = len(haplotype_lines[line[0]][3])
			haplotype_lines[line[0]][2] = line[4]
			haplotype_lines[line[0]][1] = len(line[4])
			del haplotype_lines[line[1]]
		elif line[0] == "WT":
			haplotype_lines["369"][3] += ","+line[5]
			bc_list = haplotype_lines["369"][3].split(",")
			haplotype_lines["369"][0] = len(bc_list)
			#haplotype_lines["369"][0] = len(haplotype_lines["369"][3])
			haplotype_lines["369"][2] = line[4]
			haplotype_lines["369"][1] = len(line[4])
			del haplotype_lines[line[1]]
		else:
			num_bcs = line[2]
			allele = line[4]
			len_allele = str(len(allele))
			bcs = line[5]
			haplotype_lines[line[1]] = [num_bcs,len_allele,allele,bcs]
			num_haps += 1
fixed_haplotypes_file.close()
haplotypes_file.close()

new_haplotypes_file = open("haplotypes_filtered_merged.txt","w+")
new_haplotypes_file.write("Allele\tNum_BCs\tAllele_Len\tAllele\tBCs")
num_bcs = 0
for item in haplotype_lines:
	num_bcs += int(haplotype_lines[item][0])
	haplotype_lines[item][0] = str(haplotype_lines[item][0])
	haplotype_lines[item][1] = str(haplotype_lines[item][1])
	new_haplotypes_file.write("\n"+str(item)+"\t"+"\t".join(haplotype_lines[item]))

print("Number of barcodes: " + str(num_bcs))
new_haplotypes_file.close()