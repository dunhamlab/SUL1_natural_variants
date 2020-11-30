import os
import os.path
from Bio import SeqIO
import glob
import re
import sys
#from metadata_analysis import strain_metadata

#txt_file = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-21_match_haplotypes/haplotypes_filtered_merged_sorted.txt","r")
#fasta_output_map = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-26_allele_search/191126_combined_bc_map.tsv","w+")
#fasta_output_map = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-21_np6_combined_pacbio/allele_match/191121_np6_combined_bc_map.tsv","w+")
# for line in txt_file:
# 	line = line.strip().split()
# 	line = ">"+line[0]+"\n"+line[1]+"\n"
# 	fasta_output_map.write(line)
# txt_file.close()
# fasta_output_map.close()
# fasta_output = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-21_np6_combined_pacbio/allele_match/191121_np6_combined_bc_map.tsv","r")
# fasta_output_rc = open("/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-21_np6_combined_pacbio/allele_match/191121_np6_combined_bc_map_rc.tsv","w+")
# bc_map_rc = list(SeqIO.parse(fasta_output,"fasta"))
# for item in bc_map_rc:
# 	fasta_output_rc.write(">"+str(item.id)+"\n"+str(item.seq.reverse_complement())+"\n")
# fasta_output_rc.close()
# fasta_output.close()


#sub_read_file_name = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-11-21_b_match_haplotypes/haplotypes_filtered_merged_sorted.txt"
sub_read_file_name = "../3_match_haplotypes/haplotypes_filtered_merged.txt"
# parse SUL1 allele fasta file 
ref_alleles_file = open("SUL1_ext_full.fasta", "r")
ref_alleles = list(SeqIO.parse(ref_alleles_file, "fasta"))
out_file_name = "strain_bc_map.tsv"

ref_allele_dict = {}
seq_string = ""
for item in ref_alleles:
	name = item.id
	name = name.split("_")
	if name[0] == "SACE":
		name = name[1]
	else:
		name = str(name[0])
	ref_allele_dict[name] = str(item.seq)
	seq_string = seq_string + ref_allele_dict[name]
 	#print(name + "\t" + ref_allele_dict[name])
 	#print(name)
print("Number of strains: " + str(len(ref_allele_dict.keys())))

seq_string = list(seq_string)
unique_nts = list(set(seq_string))
# ['A', 'C', 'G', 'K', 'M', 'S', 'R', 'T', 'W', 'Y']
	
# make sequence into something that can be used to grep
ref_grep_dict = {}
ambiguous_nts = ["R","Y","K","M","S","W"]
keys_w_ambiguity = []
if not os.path.isfile("SUL1_ext_full_grep.fasta.fasta"):
	ref_alleles_grep = open("SUL1_ext_full_grep.fasta","w+")
	for key in ref_allele_dict:
		grep_seq = re.sub("K","[GT]",ref_allele_dict[key])
		grep_seq = re.sub("M","[AC]", grep_seq)
		grep_seq = re.sub("S","[CG]", grep_seq)
		grep_seq = re.sub("R","[AG]", grep_seq)
		grep_seq = re.sub("W","[AT]", grep_seq)
		grep_seq = re.sub("Y","[CT]", grep_seq)
		ref_alleles_grep.write(key+"\t"+grep_seq+"\n")
		ref_grep_dict[key] = grep_seq
	ref_alleles_grep.close()
else:
	ref_alleles_grep = open("SUL1_ext_full_grep.fasta","r")
	for line in ref_alleles_grep:
		line = line.split()
		ref_grep_dict[line[0]] = line[1]
		if "[" not in line[1]:
			keys_w_ambiguity.append(line[0])
	ref_alleles_grep.close()
	print("Got grep sequences.")
	
print(len(keys_w_ambiguity))

# 	if any(nt in ref_grep_dict[key] for nt in ambiguous_nts):
# 		keys_with_ambiguity.append(key)


subassembled_reads_file = open(sub_read_file_name,"r")
sub_reads = subassembled_reads_file.readlines()
file_dict = {}
for line in sub_reads:
	if "Num_BCs" not in line:
		line2 = line.strip().split()
		file_dict[line2[0]] = line.strip() + "\t"

#strain_bc_map = {}
outfile = open(out_file_name,"w+")
outfile.write("Allele\tNum_BCs\tAllele_Len\tAllele\tBCs\tStrain\tNoAmbiguity")
for key in ref_grep_dict:
	#strain_bc_map[key] = []
	for line in sub_reads:
		if re.search(ref_grep_dict[key],line):
			line2 = line.strip().split()
			allele_index = line2[0]
			file_dict[allele_index] += key+","

for item in file_dict:
	#if any(key in file_dict[item] for key in keys_w_ambiguity) and file_dict[item][-1] == ",":
	line = file_dict[item].strip().split("\t")
	if len(line) > 5:
		line = line[5]
		if any(strain in line for strain in keys_w_ambiguity) and file_dict[item][-1] == ",":
			outfile.write("\n"+file_dict[item][:-1]+"\t"+"Yes")
		elif file_dict[item][-1] == ",":
			outfile.write("\n"+file_dict[item][:-1]+"\t"+"No")
		else:
			outfile.write("\n"+file_dict[item][:-1]+"\t"+"Error")
	else:
		outfile.write("\n"+file_dict[item][:-1]+"\t\t"+"No")

	
# read in final barcode map
# if not os.path.isfile(out_file_name):
# 	subassembled_reads_file = open(sub_read_file_name, "r")
# 	sub_reads = subassembled_reads_file.read()
# 	strain_bc_map = {}
# 	unpaired = []
# 	for key in ref_grep_dict:
# 		split_subreads = re.split(ref_grep_dict[key],sub_reads)
# 		bcs = []
# 		if len(split_subreads) > 1:
# 			for substring in split_subreads:
# 				bcs.append(substring[-11:-1])
# 		else:
# 			unpaired.append(key)
# 		strain_bc_map[key] = bcs
# 
# 	print("Number of strains without BCs: " + str(len(unpaired)))
# 
# 	subassembled_reads_file.close()
# 
# 	strain_bc_map_file = open(out_file_name, "w+")
# 	strain_bc_map_file.write("strain\tbarcodes")
# 	all_bcs = []
# 	for key in strain_bc_map:
# 		if "GGACGTTTAG" in strain_bc_map[key]:
# 			strain_bc_map[key].remove("GGACGTTTAG")
# 		all_bcs = all_bcs + strain_bc_map[key]
# 		strain_bc_map_file.write("\n"+key+"\t")
# 		strain_bc_map_file.write(",".join(strain_bc_map[key]))
# 	strain_bc_map_file.close()
# 
# 	all_bcs = list(set(all_bcs))
# 	print("Number of barcodes matched: " + str(len(all_bcs)))
# else:
# 	strain_bc_map_file = open(out_file_name, "r")
# 	strain_bc_map={}
# 	all_bcs = [] # all_bcs are actually the barcodes that matched
# 	matched = 0
# 	for line in strain_bc_map_file:
# 		line = line.split()
# 		if len(line) > 1 and "barcodes" not in line:
# 			strain = line[0]
# 			bcs = line[1].split(",")
# 			if "GGACGTTTAG" in bcs:
# 				bcs.remove("GGACGTTTAG")
# 			strain_bc_map[strain] = bcs
# 			all_bcs = all_bcs + bcs
# 			matched += 1
# 	all_bcs = list(set(all_bcs))
# 	print("Number of barcodes matched: " + str(len(all_bcs)))
# 	print("Number of strains with at least 1 bc: " + str(matched))
# 	strain_bc_map_file.close()
# 
# subassembled_reads_file = open(sub_read_file_name, "r")
# bc_read_dict = {}
# indiv_reads = list(SeqIO.parse(subassembled_reads_file, "fasta"))
# print(len(indiv_reads))
# for bc in indiv_reads:
# 	bc_read_dict[str(bc.id)] = str(bc.seq)
# 	
# unique_alleles = []
# if "barcodes" in all_bcs:
# 	all_bcs.remove("barcodes")
# for bc in all_bcs:
# 	if bc in bc_read_dict:
# 		unique_alleles.append(bc_read_dict[bc])
# 
# unique_alleles = list(set(unique_alleles))
# print("Number of unique alleles in matched barcodes: " + str(len(unique_alleles)))
# 
# # this file will include strain, number of unique alleles for this strain, and number of barcodes associated with the strain
# detailed_strain_bc_file = open("strain_alleles_details.txt","w+")
# detailed_strain_bc_file.write("Strain\tPloidy\tUnique_Alleles\tNumBCs")
# strain_bc_map_file = open(out_file_name,"r")
# for line in strain_bc_map_file:
# 	line = line.strip().split()
# 	strain = line[0]
# 	detailed_strain_bc_file.write("\n"+strain+"\t"+str(strain_metadata[strain].ploidy_num))
# 	print(strain)
# 	if len(line) > 1 and "barcodes" not in line:
# 		bcs = line[1].split(",")
# 		alleles = []
# 		for bc in bcs:
# 			alleles.append(bc_read_dict[bc])
# 		alleles = list(set(alleles))
# 		detailed_strain_bc_file.write("\t"+str(len(alleles))+"\t"+str(len(bcs))+"\t"+line[1])
# detailed_strain_bc_file.close()
# 	

# Testing for haploytpes. This is to determine the # of alleles for one strain.
# strain = sys.argv[1]
# filename = strain + "_strain_barcode_seqs.fasta"
# strain_fasta = open(filename, "w+")
# 
# strain_seqs = []
# print("Number of barcodes in " + strain + ": " + str(len(strain_bc_map[strain])))
# for bc in strain_bc_map[strain]:
# 	strain_seqs.append(bc_read_dict[bc])
# 	strain_fasta.write(">"+bc+"\n"+bc_read_dict[bc]+"\n")
# 
# filename2 = strain + "_unique_barcode_seqs.fasta"
# strain_fasta2 = open(filename2, "w+")
# strain_seqs = list(set(strain_seqs))
# i = 1
# for item in strain_seqs:
# 	strain_fasta2.write(">"+strain+"_"+str(i)+"\n"+item+"\n")
# 	i += 1
# 
# strain_fasta2.close()
# strain_fasta.close()
# subassembled_reads_file.close()
# 
# print("Number of unique sequences for " + strain + ": " + str(len(strain_seqs)))
# 


	

#for strain_name in ref_allele_dict:







	
