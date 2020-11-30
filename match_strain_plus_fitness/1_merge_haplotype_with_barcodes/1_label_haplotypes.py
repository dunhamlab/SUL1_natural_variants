# this script will help me match known haplotypes

from Bio import SeqIO

# CHANGE
bc_map_file = open("../2_allele_search/SUL1_seq_bc_rc.fasta","r")
# note: you tried this with np10 files and it turned out terrible!
#bc_map_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-21_np6_combined_pacbio/allele_match/191121_np6_combined_bc_map_rc.tsv","r")
#bc_map_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-13_match_combined_np10_alleles/191113_combined_bc_map_rc.tsv","r")
#bc_map_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-08_pb_subassembly_np3_CDS_muscle/intermediates/sub_reads_final_rc.fasta","r")

bc_map = list(SeqIO.parse(bc_map_file,"fasta"))
bc_map_file.close()

haplotypes = []
for item in bc_map:
	haplotypes.append(str(item.seq))
print(len(bc_map))
print(len(list(set(haplotypes))))

haplotype_bc_map = {}
hap_count = {}

for item in bc_map:
	bc = str(item.id)
	seq = str(item.seq)
	if seq not in haplotype_bc_map:
		haplotype_bc_map[seq] = [bc]
		hap_count[seq] = 1
	else:
		haplotype_bc_map[seq].append(bc)
		hap_count[seq] += 1

print(len(haplotype_bc_map))
count1 = 0
for item in hap_count:
	if hap_count[item] > 1:
		#print(hap_count[item])
		count1 += 1
print(count1)

haplotype_map = {}
i = 1
for item in haplotype_bc_map:
	haplotype_map[i] = item
	i += 1

haplotype_count_file = open("haplotype_counts.txt","w+")
haplotype_count_file.write("Num\tNum_Alleles\tAllele\tBCs")
haplotype_real_file =open("haplotypes_filtered.txt","w+")
haplotype_real_file.write("Num\tNum_Alleles\tAllele_Length\tAllele\tBCs")
realign_seqs_file = open("haplotypes_to_realign.txt","w+")
realign_seqs_file.write("Num\tNum_Alleles\tAllele_Length\tAllele\tBCs")
filtered_num = 1
for num in haplotype_map:
	bc_list = ",".join(haplotype_bc_map[haplotype_map[num]])
	haplotype_count_file.write("\n"+str(num)+"\t"+str(hap_count[haplotype_map[num]]) + "\t" + haplotype_map[num]+"\t"+bc_list)
	if hap_count[haplotype_map[num]] > 1:
		haplotype_real_file.write("\n"+str(filtered_num)+"\t"+str(hap_count[haplotype_map[num]]) + "\t" + str(len(haplotype_map[num])) + "\t" + haplotype_map[num]+"\t"+bc_list)
		if len(haplotype_map[num]) != 3686:
			realign_seqs_file.write("\n"+str(filtered_num)+"\t"+str(hap_count[haplotype_map[num]]) + "\t" + str(len(haplotype_map[num])) + "\t" + haplotype_map[num]+"\t"+bc_list)
		filtered_num += 1

haplotype_count_file.close()
haplotype_real_file.close()
realign_seqs_file.close()


