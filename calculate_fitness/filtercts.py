#filter out low-count bcs
import sys

sample=sys.argv[1] # sample name
read_file = open(sample+"_all_read_cts.csv", "r") # get counts file for each replicate
read_file_new = open(sample+"_cts_filtered.csv","w+") 
read_file_withbcs = open(sample+"_all_read_cts_2.csv", "w+")
for line in read_file:
	line = line.split(",")
	line2 = map(int,line[1:])
	numHighCounts = 0
	for count in line2:
		if count > 10: # if time point has more than 10 reads per time point
			numHighCounts += 1 # add one
	if numHighCounts >= 5: # if there are more than 10 reads in at least 5 time points
		if line2[0] > 10: # and time point 0 count is > 10
			read_file_new.write(",".join(line[1:])) # write file
			read_file_withbcs.write(",".join(line))
		
read_file.close()
read_file_new.close()
		