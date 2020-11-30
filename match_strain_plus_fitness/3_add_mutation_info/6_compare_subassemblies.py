# subassembly comparison between old run (1-28-20) and new run (7-30-20)

newfile = open("alleles_per_strain.txt","r")
oldfile = open("alleles_per_strain_oldrun.txt", "r")

new_strains_matched = []
for line in newfile:
	line = line.strip().split()
	new_strains_matched.append(line[0])

old_strains_matched = []
for line in oldfile:
	line = line.strip().split()
	old_strains_matched.append(line[0])

all_strains = old_strains_matched + new_strains_matched
all_strains = list(set(all_strains))

overlapping = []
in_new_not_old = []
in_old_not_new = []

for strain in all_strains:
	if strain in old_strains_matched and strain in new_strains_matched:
		overlapping.append(strain)
	elif strain in old_strains_matched and strain not in new_strains_matched:
		in_old_not_new.append(strain)
	elif strain in new_strains_matched and strain not in old_strains_matched:
		in_new_not_old.append(strain)
	else:
		print(strain)

print(len(overlapping))
print(len(in_new_not_old))
print(len(in_old_not_new))


