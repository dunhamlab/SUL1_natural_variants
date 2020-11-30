# adds strain name to fasta files

fasta_file = open("../SUL1_CDS.fasta","r")
fasta_file_fixed = open("../SUL1_ext_full_CDS.fasta","w+")

strain_names = open("strain_selected.txt","r")

names = []
for line in strain_names:
	line = line.strip()
	names.append(line)
print(len(names))

ct = 0
for line in fasta_file:
	if line[0] == ">":
		name = names[ct]
		print(ct)
		if "SACE" in name:
			#name = name.split("_")
			#name = name[1]
			fasta_file_fixed.write(">"+name+"\n")
		else:
			fasta_file_fixed.write(">"+name+"\n")
		ct += 1
	else:
		fasta_file_fixed.write(line)

print(ct)
strain_names.close()
fasta_file_fixed.close()
fasta_file.close()

