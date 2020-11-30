# this script looks at where mutations occur in pb library

#from make_pb_fastas import barcodes
from Bio import SeqIO
import os

infile1 = open("../4_allele_search/bc_strain_allele_map_mutations_fitness_median_range.txt","r")

ref_name = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/20-05-22_mutation_distribution_ref_vs_library/SUL1_ref_plusnoncoding.fasta"
ref_SUL1_fasta = open(ref_name,"r")
ref_SUL1_fasta.readline()
refseq = ref_SUL1_fasta.readline().strip()
ref_SUL1_fasta.close()

# variation_dict = {}
# for i in range(len(refseq)):
# 	variation_dict[i] = 0
ite = 0
barcodes = []
outfile = open("strain_allele_map_corrected.txt","w+")
for line in infile1:
	line1 = line.strip().split()
	bcs = line1[4]
	if "NA" not in bcs and "BCs" not in bcs:
		bcs = bcs.split(",")
		bcs = bcs[0]
		barcodes.append(bcs)
		al_file = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/20-05-22_mutation_distribution_ref_vs_library/all_pb_aligned/" + bcs + ".fasta"
		al_file = open(al_file,"r")
		al_file = SeqIO.parse(al_file,"fasta")
		for item in al_file:
			if "SUL1_ref" in str(item.id):
				refalignment = str(item.seq)
			else:
				queryalignment = str(item.seq)
		counter = 0
		num_changes = 0
		mutations = []
		for i in range(len(refalignment)):
			if refalignment[i] != queryalignment[i]:
				num_changes += 1
				#mut = refalignment[i] + str(i+counter) + queryalignment[i]
				#mutations.append(mut)
				if refalignment[i] == "-":
					mut = "ins"+str(i+counter-843)+queryalignment[i]
					counter -= 1
				elif queryalignment[i] == "-":
					mut = "del"+str(i+counter-843)+refalignment[i]
				else:
					mut = refalignment[i] + str(i+counter-843) + queryalignment[i]
				mutations.append(mut)
		mutations = mutations[:-1]
		if len(mutations) > 0:
			line1[7] = ",".join(mutations)
			line1[8] = str(num_changes)
		else:
			line1[7] = "NA"
			line1[8] = "0"
		line2 = "\t".join(line1)+"\n"
		outfile.write(line2)
	else:
		if "BCs" in line:
			line1 = line.strip().split()
			line1[8] = "NumNtChanges"
			line = "\t".join(line1)
			outfile.write(line+"\n")
		elif "NA" in line:
			line1 = line.strip().split()
			new_fasta_name = "pb_fasta/allele" + line[0] + ".fasta"
			new_fasta = open(new_fasta_name,"w+")
			new_fasta.write(">allele" + line1[0] + "\n" + line1[3]+"\n")
			new_fasta.close()
			outfile_name = "pb_aligned/allele" + line[0] + ".fasta"
			#cmd = "needle " + new_fasta_name + " " + ref_name + " -gapopen 10 -gapextend 0.5 -outfile " + outfile_name + " -aformat fasta"
			#os.system(cmd)
			al_file = open(outfile_name,"r")
			al_file = SeqIO.parse(al_file,"fasta")
			for item in al_file:
				if "SUL1_ref" in str(item.id):
					refalignment = str(item.seq)
				else:
					queryalignment = str(item.seq)
			counter = 0
			num_changes = 0
			mutations = []
			for i in range(len(refalignment)):
				if refalignment[i] != queryalignment[i]:
					num_changes += 1
					#mut = refalignment[i] + str(i+counter) + queryalignment[i]
					#mutations.append(mut)
					if refalignment[i] == "-":
						mut = "ins"+str(i+counter-843)+queryalignment[i]
						counter -= 1
					elif queryalignment[i] == "-":
						mut = "del"+str(i+counter-843)+refalignment[i]
					else:
						mut = refalignment[i] + str(i+counter-843) + queryalignment[i]
					mutations.append(mut)
			al_file.close()
			#new_outfile = open
			line1[7] = ",".join(mutations)
			line1[8] = str(num_changes)
			line2 = "\t".join(line1)+"\n"
			outfile.write(line2)
			ite += 1
			print(ite)
		else:
			outfile.write(line)
	
	
infile1.close()




# ite = 0
# for bc in barcodes:
# 	al_file = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/20-05-22_mutation_distribution_ref_vs_library/all_pb_aligned/" + bc + ".fasta"
# 	al_file = open(al_file,"r")
# 	al_file = SeqIO.parse(al_file,"fasta")
# 	for item in al_file:
# 		if "SUL1_ref" in str(item.id):
# 			refalignment = str(item.seq)
# 		else:
# 			queryalignment = str(item.seq)
# 	counter=0
# 	num_changes = 0
# 	for i in range(len(refalignment)):
# 		if refalignment[i] != queryalignment[i]:
# 			num_changes += 1
# 			variation_dict[i+counter] += 1
# 			if refalignment[i] == "-":
# 				counter -= 1
# 	if ite % 200 == 0 :
# 		print(ite)
# 	ite += 1
# 	al_file.close()
# pb_var_file = open("pbs_mutation_distribution_updated_withnoncoding.txt","w+")
# for i in range(len(refseq)):
# 	pb_var_file.write(str(i-843) + "\t" + str(variation_dict[i])+"\n")
# pb_var_file.close()






