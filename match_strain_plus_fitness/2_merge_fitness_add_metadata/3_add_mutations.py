# this python script adds annotation for which mutations occurred. this is is done on 191126_updated_strain_bc_map

from Bio.Seq import Seq
in_file = open("updated_strain_bc_map.tsv","r")
out_file = open("bc_strain_allele_map_plusmutations.txt","w+")

ref_allele = "ATGTCACGTAAGAGCTCGACTGAATATGTGCATAATCAGGAGGATGCTGATATCGAAGTATTTGAATCAGAATACCGCACATATAGGGAATCTGAGGCGGCAGAAAACAGAGACGGACTTCACAATGGTGATGAGGAAAATTGGAAGGTTAATAGTAGTAAGCAGAAATTTGGGGTAACGAAAAATGAGCTATCAGATGTCCTGTACGATTCCATTCCAGCGTATGAAGAGAGCACAGTCACTTTGAAGGAGTACTATGATCATTCTATCAAAAACAATCTAACTGCGAAATCGGCAGGAAGTTACCTCGTATCTCTTTTTCCTATTATAAAATGGTTTCCTCATTATAACTTTACGTGGGGCTATGCTGATTTAGTGGCAGGAATTACAGTTGGCTGCGTACTCGTGCCCCAATCTATGTCATACGCACAAATCGCTAGTTTATCTCCTGAATATGGTTTGTATTCCTCCTTTATTGGTGCGTTTATATATTCTTTGTTTGCCACATCGAAAGATGTTTGTATTGGTCCGGTCGCTGTAATGTCACTACAAACTGCCAAAGTCATTGCTGAAGTTCTAAAAAAATATCCCGAAGACCAGACAGAAGTTACAGCTCCTATCATTGCAACTACCCTTTGTTTGCTTTGTGGGATTGTCGCCACTGGGTTGGGTATACTGCGTTTAGGCTTTTTAGTGGAACTTATTTCTCTAAATGCTGTTGCTGGCTTCATGACCGGTTCCGCATTTAACATCATCTGGGGTCAAATTCCGGCTCTCATGGGATACAACTCATTAGTGAATACCAGAGAAGCAACGTATAAGGTTGTAATTAACACTCTGAAACATTTACCAAACACAAAGTTAGACGCCGTTTTTGGCTTGATTCCGTTGGTAATCCTCTATGTATGGAAATGGTGGTGTGGTACATTTGGTATAACTTTGGCAGATAGATATTATCGAAATCAACCAAAGGTAGCAAATAGACTGAAATCCTTCTATTTCTATGCACAAGCTATGAGAAATGCCGTCGTCATAGTAGTTTTTACTGCCATATCGTGGAGCATAACAAGAAACAAATCTTCAAAAGACCGTCCAATCAGTATTCTGGGTACAGTTCCCTCGGGCTTAAATGAGGTGGGAGTTATGAAAATCCCAGACGGTCTGCTATCTAATATGAGTTCAGAAATACCTGCTTCAATTATCGTTCTGGTGTTAGAACACATCGCTATTTCAAAATCCTTTGGTAGAATTAACGACTACAAGGTTGTCCCTGACCAAGAACTTATTGCGATTGGTGTGACAAATTTGATAGGGACATTTTTTCACTCATATCCAGCAACTGGGTCATTTTCCAGATCTGCTTTGAAAGCAAAATGTAACGTGCGCACTCCGTTTTCTGGGGTATTCACTGGCGGTTGCGTTCTATTAGCACTTTATTGTTTAACTGACGCCTTCTTTTTCATTCCTAAAGCGACACTATCGGCGGTTATTATTCATGCTGTTTCTGATTTGCTGACTTCTTACAAAACCACCTGGACCTTCTGGAAGACCAACCCGTTAGATTGTATCTCATTTATCGTTACAGTGTTCATCACAGTATTTTCATCCATTGAAAATGGTATATATTTTGCAATGTGTTGGTCATGTGCAATGTTACTATTGAAACAGGCTTTCCCTGCTGGTAAATTCCTTGGTCGTGTTGAGGTGGCAGAAGTATTGAACCCAACAGTACAAGAGGATATTGATGCTGTGATATCATCTAATGAATTACCTAATGAACTGAATAAACAGGTTAAGTCTACTGTTGAGGTTTTACCAGCCCCAGAGTATAAGTTTAGCGTAAAGTGGGTTCCGTTCGATCATGGATACTCAAGAGAATTGAATATCAATACCACAGTTCGGCCTCCTCCACCAGGTGTCATAGTCTATCGTTTGGGTGATAGCTTTACTTACGTGAACTGCTCAAGGCATTATGACATTATATTTGATCGTATTAAGGAAGAAACAAGGCGAGGCCAACTTATAACCTTAAGGAAAAAGTCAGACCGTCCATGGAATGATCCTGGTGAATGGAAAATGCCAGATTCTTTGAAATCACTATTTAAATTTAAACGTCATTCAGCAACAACGAATAGTGACCTACCGATATCGAATGGAAGCAGTAACGGAGAAACATATGAAAAGCCGCTACTGAAAGTCGTCTGCCTGGATTTTTCCCAAGTTGCTCAAGTGGATTCAACCGCTGTTCAAAGCCTGGTTGATCTGAGAAAAGCTGTGAATAGGTATGCGGATAGACAAGTCGAATTCCATTTTGCCGGAATTATATCTCCATGGATCAAAAGAAGTCTTTTGAGTGTTAAATTCGGAACTACAAATGAGGAATATAGTGACGACTCTATTATCGCTGGCCATTCTAGTTTTCACGTTGCAAAAGTTTTGAAGGATGATGTGGATTATACTGATGAAGACAGCCGTATAAGCACATCTTACAGTAACTATGAAACATTATGTGCTGCAACTGGGACAAATTTACCGTTTTTTCATATCGATATACCCGATTTTTCTAAATGGGACGTTTAG"
print("length of ref: " + str(len(ref_allele)))
ref_aa = Seq(ref_allele)
ref_aa = ref_aa.translate()


for line in in_file:
	line = line.strip().split()
	mutations=[]
	allele = line[3]
	numIdentical = 0
	aa_mutations=[]
	aaIdentical = 0
	if line[0] != "Allele":
		print("length of allele: " + str(len(allele)))
		allele_aa = Seq(allele).translate()
		for i in range(len(ref_allele)):
			if ref_allele[i] == allele[i]:
				numIdentical += 1
			else:
				mut = ref_allele[i]+str(i+1)+allele[i]
				mutations.append(mut)
		#perIdentity = numIdentical / len(allele)
		for i in range(len(ref_aa)):
			if ref_aa[i] == allele_aa[i]:
				aaIdentical += 1
			else:
				mut = ref_aa[i]+str(i+1)+allele_aa[i]
				aa_mutations.append(mut)
		if len(line) == 7:
			out_file.write(("\t".join(line)))
			if len(aa_mutations) > 0:
				out_file.write("\t"+",".join(mutations)+"\t"+str(numIdentical)+"\t"+",".join(aa_mutations)+"\t"+str(aaIdentical)+"\n")
			else:
				out_file.write("\t"+",".join(mutations)+"\t"+str(numIdentical)+"\t"+"NA"+"\t"+str(aaIdentical)+"\n")
		else:
			out_file.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\tNA\t"+line[4]+"\t"+line[5]+"\t"+",".join(mutations)+"\t"+str(numIdentical)+"\t"+",".join(aa_mutations)+"\t"+str(aaIdentical)+"\n")
	else:
		out_file.write("\t".join(line)+"\tNtMutations\tNumIdentical\tAAMutations\tNumAAIdentical"+"\n")





# for line in in_file:
# 	line = line.strip().split()
# 	mutations=[]
# 	allele = line[3]
# 	numIdentical = 0
# 	allele_aa = Seq(allele).translate()
# 	aa_mutations=[]
# 	aaIdentical = 0
# 	if line[0] != "Allele":
# 		print("length of allele: " + str(len(allele)))
# 		for i in range(len(ref_allele)):
# 			if ref_allele[i] == allele[i]:
# 				numIdentical += 1
# 			else:
# 				mut = ref_allele[i]+str(i)+allele[i]
# 				mutations.append(mut)
# 		#perIdentity = numIdentical / len(allele)
# 		for i in range(len(ref_aa)):
# 			if ref_aa[i] == allele_aa[i]:
# 				aaIdentical += 1
# 			else:
# 				mut = ref_aa[i]+str(i)+allele_aa[i]
# 				aa_mutations.append(mut)
# 		if len(line) == 7:
# 			out_file.write(("\t".join(line)))
# 			out_file.write("\t"+",".join(mutations)+"\t"+str(numIdentical)+"\t"+",".join(aa_mutations)+"\t"+str(aaIdentical)+"\n")
# 		else:
# 			out_file.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\tNA\t"+line[4]+"\t"+line[5]+"\t"+",".join(mutations)+"\t"+str(numIdentical)+"\t"+",".join(aa_mutations)+"\t"+str(aaIdentical)+"\n")
# 	else:
# 		out_file.write("\t".join(line)+"\tNtMutations\tNumIdentical\tAAMutations\tNumAAIdentical"+"\n")

in_file.close()
out_file.close()

#out_file = open("200109_bc_strain_allele_map_plusmutations.txt","a")
