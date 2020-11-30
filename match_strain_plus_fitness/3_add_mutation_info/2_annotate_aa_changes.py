# this script adds amino acid changes to file

from Bio import SeqIO
from Bio.Seq import Seq
import re
import os

refname = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/20-05-22_mutation_distribution_ref_vs_library/SUL1_ref_plusnoncoding.fasta"
reffile = open(refname,"r")
reffile.readline()
refseq = reffile.readline()
print(refseq)
CDS = refseq[844:2580+844]
Seq_CDS = Seq(CDS)
ref_aa = str(Seq_CDS.translate())
ref_aa_file_name = "ref_SUL1_aa.fasta"
ref_aa_file = open("ref_SUL1_aa.fasta","w+")
ref_aa_file.write(">SUL1_aa_ref\n"+ref_aa+"\n")
ref_aa_file.close()
# refseq_dict = {}
# for i in range(len(refseq)):
# 	refseq_dict[i-843] = refseq[i]

infile = open("strain_allele_map_corrected.txt","r")
outfile = open("strain_allele_map_corrected_2.txt","w+")
ite = 1
for line in infile:
	print(ite)
	ite += 1
	line = line.strip().split()
	if "Allele" not in line:
		temprefseq = CDS
		mutations = line[7].split(",")
		if "NA" not in mutations:
			counter = 0
			for mut in mutations:
				print(mut)
				if "ins" in mut:
					#mut1 = mut[0:3]
					ins_location = int(mut[3:-1])
					if ins_location >= 0 and ins_location < 2580:
						added_nt = mut[-1]
						seq1 = temprefseq[0:ins_location+counter]
						seq2 = added_nt
						seq3 = temprefseq[ins_location+counter:]
						temprefseq = seq1+seq2+seq3
						counter += 1
				elif "del" in mut:
					#mut1 = [0:3]
					del_location = int(mut[3:-1])
					if del_location >= 0 and del_location < 2580:
						seq1 = temprefseq[0:del_location+counter]
						seq2 = temprefseq[del_location+counter+1:]
						temprefseq = seq1+seq2
						counter -= 1
				else:
					mut_location = int(mut[1:-1])-1
					print(mut_location)
					if mut_location >= 0 and mut_location < 2580:
						changed_nt = mut[-1]
						seq1 = temprefseq[0:mut_location+counter]
						seq2 = changed_nt
						seq3 = temprefseq[mut_location+1+counter:]
						temprefseq = seq1+seq2+seq3
			#print(temprefseq)
			temprefseq = Seq(temprefseq)
			query_aa = str(temprefseq.translate())
			query_aa_file_name = "pb_aa_fasta/allele"+line[0]+".fasta"
			#query_aa_file = open(query_aa_file_name,"w+")
			#query_aa_file.write(">allele"+line[0]+"\n"+str(query_aa)+"\n")
			#query_aa_file.close()
			outfile_name = "pb_aa_aligned/allele"+line[0]+".fasta"
			#cmd = "needle " + query_aa_file_name + " " + ref_aa_file_name + " -gapopen 10 -gapextend 0.5 -outfile " + outfile_name + " -aformat fasta"
			#os.system(cmd)
			
			#alignment_file = open(outfile_name,"r") # original
			#alignment = SeqIO.parse(alignment_file,"fasta")
			#for item in alignment:
			#	if str(item.id) == "SUL1_aa_ref":
			#		refseq = str(item.seq)
			#	else:
			#		queryseq = str(item.seq)
			
			all_aa_changes = []
			num_changes = 0
			print(query_aa)
			if len(query_aa) >= len(ref_aa):
				for i in range(len(ref_aa)):
					if ref_aa[i] != query_aa[i]:
						aa_change = ref_aa[i] + str(i+1) + query_aa[i]
						all_aa_changes.append(aa_change)
						num_changes += 1
			elif len(query_aa) < len(ref_aa):
				for i in range(len(query_aa)):
					if ref_aa[i] != query_aa[i]:
						aa_change = ref_aa[i] + str(i+1) + query_aa[i]
						all_aa_changes.append(aa_change)
						num_changes += 1
			if len(all_aa_changes) > 0:
				all_aa_changes = ",".join(all_aa_changes)
			else:
				all_aa_changes = "NA"
							
# 			all_aa_changes = []
# 			num_changes = 0
# 			for i in range(len(refseq)):
# 				if refseq[i] != queryseq[i]:
# 					aa_change = refseq[i] + str(i+1) + queryseq[i]
# 					all_aa_changes.append(aa_change)
# 					num_changes += 1
# 			if len(all_aa_changes) > 0:
# 				all_aa_changes = ",".join(all_aa_changes)
# 			else:
# 				all_aa_changes = "NA"

			line[9] = all_aa_changes
			line[10] = str(num_changes)
			line.append(str(temprefseq))
		else:
			line[9] = "NA"
			line[10] = "0"
			line.append("NA")
	else:
		line[10] = "NumAAChanges"
	outfile.write("\t".join(line)+"CDS\n")
			
			
# need to manually check!!
			
			
			
			