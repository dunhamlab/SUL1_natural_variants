# this script analyzes the metadata for the matched alleles


# with open("strain_metadata.txt", encoding="utf-8") as f:
#     print(f.read())
    
# with open("strain_metadata.txt", "rb") as file:
#      data = file.read(8)
#      data = str(data)
#      print(data)
# with open("strain_metadata_2.txt", "w+") as f:
#    f.write(" ".join(map(str,data)))
#    f.write("\n")
# 
metadata_file = open("strain_metadata_2.txt","r") # this is the supplemental file from the 1,011 genomes paper with all the metadata for each strain


class Strain:
	
	name = "NA"
	def __init__(self,name):
		self.name = name
	
	def isolate_name(self,isolate):
		self.isolate = isolate
	
	def isolated_from(self,place):
		self.isolation_place = place
	
	def eco_origin(self,eco):
		self.eco = eco
	
	def geo_origin(self,geo):
		self.geo = geo
		
	def HO_del(self,HO):
		self.HO = HO
		
	def plasmid_CN(self,plasmidCN):
		self.plasmidCN = plasmidCN
		
	def ploidy_num(self,ploidy):
		self.ploidy = ploidy
		
	aneuploidy = False
	def change_aneuploid(self,new_aneuploidy):
		self.aneuploidy = new_aneuploidy
		
	def zygosity(self,zygosity):
		self.zygosity = zygosity
		
	def coverage(self,coverage):
		self.coverage = coverage
		
	def clade_name(self,clade):
		self.clade = clade
		
metadata_file = metadata_file.readlines()
metadata_file = metadata_file[1:]

strain_metadata = {}
for line in metadata_file:
	line = str(line)
	strain = line.strip().split("\t")
	print(line)
	if len(strain[1]) > 3:
		strain[1] = strain[1][5:8]
	myStrain = Strain(strain[1])
	myStrain.isolate_name = strain[0]
	myStrain.isolated_from = strain[2]
	myStrain.eco_origin = strain[3]
	myStrain.geo_origin = strain[4]
	if strain[8] == "Xn":
		myStrain.ploidy_num = None
	else:
		myStrain.ploidy_num = int(strain[8])
	if strain[9] != "euploid":
		myStrain.change_aneuploid(True)
	myStrain.zygosity = strain[10]
	myStrain.coverage = strain[14]
	myStrain.clade_name = strain[15]
	strain_metadata[myStrain.name] = myStrain
	
for item in strain_metadata:
	print("Name: " + strain_metadata[item].name)
	print("Isolate Name: " + strain_metadata[item].isolate_name)
	print("Isolated From: " + strain_metadata[item].isolated_from)
	print("Eco: " + strain_metadata[item].eco_origin)
	print("Geo: " + strain_metadata[item].geo_origin)
	print("Ploidy: " + str(strain_metadata[item].ploidy_num))
	print("Aneuploid? " + str(strain_metadata[item].aneuploidy))
	print("Zygosity: " + strain_metadata[item].zygosity)
	print("Coverage: " + strain_metadata[item].coverage)
	print("Clade: " + strain_metadata[item].clade_name)
	print("\n")
	
#print(type(strain_metadata["CEM"].coverage))
#print(strain_metadata["CEM"].coverage*2)
	
	
	