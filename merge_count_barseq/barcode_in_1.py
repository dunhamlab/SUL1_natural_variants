import math
#import numpy as np
import random

# these functions are for dealing with very small values
# redo of barseq analysis

def eexp(x):
    if x == float('nan'):
        return 0
    else:
        return math.exp(x)


def eln(x):
    if x == 0:
        return float('nan')
    elif x > 0:
        return math.log(x)
    else:
        print "error"


def elnsum(x, y):  # x and y should already be e^x or e^y (in exponent form)
    if x == float('nan') or y == float('nan'):
        if x == math.nan:
            return y
        else:
            return x
    else:
        if x > y:
            return x + eln(1+math.exp(y-x))
        else:
            return y + eln(1+math.exp(x-y))


def bc_counts(fileName, outDict):  # get number of reads for each barcode, including ones with only 1 read
    line = fileName.readline()
    # outDict = {}
    while line:
        if line[0] == "@":
            line2 = fileName.readline()
            if not line2[0:10] in outDict:
                outDict[line2[0:10]] = 1
            else:
                outDict[line2[0:10]] += 1
            fileName.readline()
            fileName.readline()
            line = fileName.readline()
    return outDict


def unique_bcs(myDict):  # print number of unique barcodes over the total number of reads
    i = 0
    for key in myDict:
        if myDict[key] > 1:
            i += 1
    print "unique barcodes: ", i
    print "number of reads: ", sum(myDict.values())


def normalize(txDict, readCt, normDict, templist):
    print "unique bcs, including freq = 1: ", len(txDict)  # number of unique reads, including those with only 1 bc
    for key in list(txDict.keys()):  # get rid of bc's with only 1 read
        if txDict[key] < 10: #== 1:
            del txDict[key]
    for key in txDict:
        templist.append(txDict[key])
    print "smallest number of reads: ", min(templist)
    print "biggest number of reads: ", max(templist)
    print "average number of reads: ", sum(templist)/len(templist)
    print "unique bcs, excluding freq < 10: ", len(txDict)  # number of unique bcs >1 read
    for key in list(txDict.keys()):
        if key in normDict:
            txDict[key] = eln(txDict[key]) - eln(readCt)
            txDict[key] = txDict[key] - normDict[key]
        else:
            del txDict[key]
    print "unique bcs >10 read, also detected in T0: ", len(txDict)  # number of unique bcs >1 read, without ones that weren't detected in T0


def writehist(myDict, file):
    for key in myDict:
        file.write(str(key))
        file.write("\t")
        file.write(str(myDict[key]))
        file.write("\n")
    print "file written"


def writePoints(myDict, myFile):
    for key in myDict:
        myFile.write(key)
        i = 0
        for i in range(len(myDict[key])):
            myFile.write("\t")
            myFile.write(str(myDict[key][i]))
        myFile.write("\n")


sample = "61" # or change to any sample name
readFolder = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-04-30_barseq_fitness/paired_fastqs/" # these were generated from merging/counting barseq files
cwd = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-06-03_barseq_reanalysis/"
outputLoc = "/net/dunham/vol2/Cindy/SUL1_alleleLibraryCompetition/results/19-06-03_barseq_reanalysis/output/"

# these are values that are to be used as normalization
print "T0/setting normalization stuff"
#T0 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T0_S3_R1_001.fastq", "r")
timepoint = "T0"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T0 = open(sampleFileName,"r")
T0_dict = {}
bc_counts(T0, T0_dict)  # count # of reads per barcode
T0_readct = {}
for key in T0_dict:
	T0_readct[key] = T0_dict[key]
unique_bcs(T0_dict)  # get an idea of # of unique BCs and total # of reads
T0.close()
T0_reads = sum(T0_dict.values())  # total number of reads
print "unique bcs, including freq = 1: ", len(T0_dict)  # number of unique barcodes total, including ones with 1 read
for key in list(T0_dict.keys()):  # remove BCs that only have 1 read. probably sequencing error or noise.
    if T0_dict[key] == 1:
        del T0_dict[key]
readsHist_T0 = []
for key in list(T0_dict.keys()):
    readsHist_T0.append(T0_dict[key])
print "unique bcs, excluding freq <= 1: ", len(T0_dict)  # number of unique barcodes total, excluding ones with 1 read
for key in T0_dict:  # normalize to number of reads
    T0_dict[key] = eln(T0_dict[key]) - eln(T0_reads)
#  write log values to file
T0_hist_out = open(outputLoc+"T0_hist.txt", "w")
writehist(T0_dict, T0_hist_out)
T0_hist_out.close()
T0_rawReads = open(outputLoc+"T0_rawReads.txt", "w")
for num in readsHist_T0:
    T0_rawReads.write(str(num))
    T0_rawReads.write("\n")
T0_rawReads.close()
print "make sure nothing changed: ", len(T0_dict), "\n"


print "T1"
timepoint = "T1"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T1 = open(sampleFileName,"r")
#T1 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T1_S4_R1_001.fastq", "r")
T1_dict = {}
bc_counts(T1, T1_dict)  # count # of reads per barcode
T1_readct = {}
for key in T1_dict:
	T1_readct[key] = T1_dict[key]
unique_bcs(T1_dict)  # get an idea of # of unique BCs and total # of reads
T1.close()
T1_reads = sum(T1_dict.values())  # total number of reads
readsHist_T1 = []
normalize(T1_dict, T1_reads, T0_dict, readsHist_T1)
#normalize(T1_dict, T1_reads, T0_dict, readsHist_T1, T1_dict_filtered)  # all values should be normalized now?
# write values to file
T1_hist_out = open(outputLoc+"T1_hist.txt", "w")
writehist(T1_dict, T1_hist_out) # this is for normalized reads
# writehist(T1_dict_filtered, T1_hist_out) # this is for raw read counts
T1_hist_out.close()
T1_rawReads = open(outputLoc+"T1_rawReads.txt", "w")
for num in readsHist_T1:
    T1_rawReads.write(str(num))
    T1_rawReads.write("\n")
T1_rawReads.close()
print "make sure nothing changed: ", len(T1_dict), "\n"


print "T2"
#T2 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T2_S5_R1_001.fastq", "r")
timepoint = "T2"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T2 = open(sampleFileName,"r")
T2_dict = {}
bc_counts(T2, T2_dict)  # count # of reads per barcode
T2_readct = {}
for key in T2_dict:
	T2_readct[key] = T2_dict[key]
unique_bcs(T2_dict)  # get an idea of # of unique BCs and total # of reads
T2.close()
T2_reads = sum(T2_dict.values())  # total number of reads
readsHist_T2 = []
T2_dict_filtered = {}
normalize(T2_dict, T2_reads, T0_dict, readsHist_T2)  # all values should be normalized now?
# write values to file
T2_hist_out = open(outputLoc+"T2_hist.txt", "w")
writehist(T2_dict, T2_hist_out) # this is for normalized reads
# writehist(T2_dict_filtered, T2_hist_out) # this is for raw read counts
T2_hist_out.close()
T2_rawReads = open(outputLoc+"T2_rawReads.txt", "w")
for num in readsHist_T2:
    T2_rawReads.write(str(num))
    T2_rawReads.write("\n")
T2_rawReads.close()
print "make sure nothing changed: ", len(T2_dict), "\n"

print "T3"
#T3 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T3_S5_R1_001.fastq", "r")
timepoint = "T3"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T3 = open(sampleFileName,"r")
T3_dict = {}
bc_counts(T3, T3_dict)  # count # of reads per barcode
T3_readct = {}
for key in T3_dict:
	T3_readct[key] = T3_dict[key]
unique_bcs(T3_dict)  # get an idea of # of unique BCs and total # of reads
T3.close()
T3_reads = sum(T3_dict.values())  # total number of reads
readsHist_T3 = []
T3_dict_filtered = {}
normalize(T3_dict, T3_reads, T0_dict, readsHist_T3)  # all values should be normalized now?
# write values to file
T3_hist_out = open(outputLoc+"T3_hist.txt", "w")
writehist(T3_dict, T3_hist_out) # this is for normalized reads
# writehist(T3_dict_filtered, T3_hist_out) # this is for raw read counts
T3_hist_out.close()
T3_rawReads = open(outputLoc+"T3_rawReads.txt", "w")
for num in readsHist_T3:
    T3_rawReads.write(str(num))
    T3_rawReads.write("\n")
T3_rawReads.close()
print "make sure nothing changed: ", len(T3_dict), "\n"


print "T4"
#T4 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T4_S5_R1_001.fastq", "r")
timepoint = "T4"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T4 = open(sampleFileName,"r")
T4_dict = {}
bc_counts(T4, T4_dict)  # count # of reads per barcode
T4_readct = {}
for key in T4_dict:
	T4_readct[key] = T4_dict[key]
unique_bcs(T4_dict)  # get an idea of # of unique BCs and total # of reads
T4.close()
T4_reads = sum(T4_dict.values())  # total number of reads
readsHist_T4 = []
T4_dict_filtered = {}
normalize(T4_dict, T4_reads, T0_dict, readsHist_T4)  # all values should be normalized now?
# write values to file
T4_hist_out = open(outputLoc+"T4_hist.txt", "w")
writehist(T4_dict, T4_hist_out) # this is for normalized reads
# writehist(T4_dict_filtered, T4_hist_out) # this is for raw read counts
T4_hist_out.close()
T4_rawReads = open(outputLoc+"T4_rawReads.txt", "w")
for num in readsHist_T4:
    T4_rawReads.write(str(num))
    T4_rawReads.write("\n")
T4_rawReads.close()
print "make sure nothing changed: ", len(T4_dict), "\n"

print "T5"
#T5 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T5_S5_R1_001.fastq", "r")
timepoint = "T5"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T5 = open(sampleFileName,"r")
T5_dict = {}
bc_counts(T5, T5_dict)  # count # of reads per barcode
T5_readct = {}
for key in T5_dict:
	T5_readct[key] = T5_dict[key]
unique_bcs(T5_dict)  # get an idea of # of unique BCs and total # of reads
T5.close()
T5_reads = sum(T5_dict.values())  # total number of reads
readsHist_T5 = []
T5_dict_filtered = {}
normalize(T5_dict, T5_reads, T0_dict, readsHist_T5)  # all values should be normalized now?
# write values to file
T5_hist_out = open(outputLoc+"T5_hist.txt", "w")
writehist(T5_dict, T5_hist_out) # this is for normalized reads
# writehist(T5_dict_filtered, T5_hist_out) # this is for raw read counts
T5_hist_out.close()
T5_rawReads = open(outputLoc+"T5_rawReads.txt", "w")
for num in readsHist_T5:
    T5_rawReads.write(str(num))
    T5_rawReads.write("\n")
T5_rawReads.close()
print "make sure nothing changed: ", len(T5_dict), "\n"

print "T6"
#T6 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T6_S5_R1_001.fastq", "r")
timepoint = "T6"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T6 = open(sampleFileName,"r")
T6_dict = {}
bc_counts(T6, T6_dict)  # count # of reads per barcode
T6_readct = {}
for key in T6_dict:
	T6_readct[key] = T6_dict[key]
unique_bcs(T6_dict)  # get an idea of # of unique BCs and total # of reads
T6.close()
T6_reads = sum(T6_dict.values())  # total number of reads
readsHist_T6 = []
T6_dict_filtered = {}
normalize(T6_dict, T6_reads, T0_dict, readsHist_T6)  # all values should be normalized now?
# write values to file
T6_hist_out = open(outputLoc+"T6_hist.txt", "w")
writehist(T6_dict, T6_hist_out) # this is for normalized reads
# writehist(T6_dict_filtered, T6_hist_out) # this is for raw read counts
T6_hist_out.close()
T6_rawReads = open(outputLoc+"T6_rawReads.txt", "w")
for num in readsHist_T6:
    T6_rawReads.write(str(num))
    T6_rawReads.write("\n")
T6_rawReads.close()
print "make sure nothing changed: ", len(T6_dict), "\n"


print "T7"
#T7 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T7_S5_R1_001.fastq", "r")
timepoint = "T7"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T7 = open(sampleFileName,"r")
T7_dict = {}
bc_counts(T7, T7_dict)  # count # of reads per barcode
T7_readct = {}
for key in T7_dict:
	T7_readct[key] = T7_dict[key]
unique_bcs(T7_dict)  # get an idea of # of unique BCs and total # of reads
T7.close()
T7_reads = sum(T7_dict.values())  # total number of reads
readsHist_T7 = []
T7_dict_filtered = {}
normalize(T7_dict, T7_reads, T0_dict, readsHist_T7)  # all values should be normalized now?
# write values to file
T7_hist_out = open(outputLoc+"T7_hist.txt", "w")
writehist(T7_dict, T7_hist_out) # this is for normalized reads
# writehist(T7_dict_filtered, T7_hist_out) # this is for raw read counts
T7_hist_out.close()
T7_rawReads = open(outputLoc+"T7_rawReads.txt", "w")
for num in readsHist_T7:
    T7_rawReads.write(str(num))
    T7_rawReads.write("\n")
T7_rawReads.close()
print "make sure nothing changed: ", len(T7_dict), "\n"


print "T8"
#T8 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T8_S5_R1_001.fastq", "r")
timepoint = "T8"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T8 = open(sampleFileName,"r")
T8_dict = {}
bc_counts(T8, T8_dict)  # count # of reads per barcode
T8_readct = {}
for key in T8_dict:
	T8_readct[key] = T8_dict[key]
unique_bcs(T8_dict)  # get an idea of # of unique BCs and total # of reads
T8.close()
T8_reads = sum(T8_dict.values())  # total number of reads
readsHist_T8 = []
T8_dict_filtered = {}
normalize(T8_dict, T8_reads, T0_dict, readsHist_T8)  # all values should be normalized now?
# write values to file
T8_hist_out = open(outputLoc+"T8_hist.txt", "w")
writehist(T8_dict, T8_hist_out) # this is for normalized reads
# writehist(T8_dict_filtered, T8_hist_out) # this is for raw read counts
T8_hist_out.close()
T8_rawReads = open(outputLoc+"T8_rawReads.txt", "w")
for num in readsHist_T8:
    T8_rawReads.write(str(num))
    T8_rawReads.write("\n")
T8_rawReads.close()
print "make sure nothing changed: ", len(T8_dict), "\n"


print "T9"
#T9 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T9_S5_R1_001.fastq", "r")
timepoint = "T9"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T9 = open(sampleFileName,"r")
T9_dict = {}
bc_counts(T9, T9_dict)  # count # of reads per barcode
T9_readct = {}
for key in T9_dict:
	T9_readct[key] = T9_dict[key]
unique_bcs(T9_dict)  # get an idea of # of unique BCs and total # of reads
T9.close()
T9_reads = sum(T9_dict.values())  # total number of reads
readsHist_T9 = []
T9_dict_filtered = {}
normalize(T9_dict, T9_reads, T0_dict, readsHist_T9)  # all values should be normalized now?
# write values to file
T9_hist_out = open(outputLoc+"T9_hist.txt", "w")
writehist(T9_dict, T9_hist_out) # this is for normalized reads
# writehist(T9_dict_filtered, T9_hist_out) # this is for raw read counts
T9_hist_out.close()
T9_rawReads = open(outputLoc+"T9_rawReads.txt", "w")
for num in readsHist_T9:
    T9_rawReads.write(str(num))
    T9_rawReads.write("\n")
T9_rawReads.close()
print "make sure nothing changed: ", len(T9_dict), "\n"

print "T10"
#T10 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T10_S5_R1_001.fastq", "r")
timepoint = "T10"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T10 = open(sampleFileName,"r")
T10_dict = {}
bc_counts(T10, T10_dict)  # count # of reads per barcode
T10_readct = {}
for key in T10_dict:
	T10_readct[key] = T10_dict[key]
unique_bcs(T10_dict)  # get an idea of # of unique BCs and total # of reads
T10.close()
T10_reads = sum(T10_dict.values())  # total number of reads
readsHist_T10 = []
T10_dict_filtered = {}
normalize(T10_dict, T10_reads, T0_dict, readsHist_T10)  # all values should be normalized now?
# write values to file
T10_hist_out = open(outputLoc+"T10_hist.txt", "w")
writehist(T10_dict, T10_hist_out) # this is for normalized reads
# writehist(T10_dict_filtered, T10_hist_out) # this is for raw read counts
T10_hist_out.close()
T10_rawReads = open(outputLoc+"T10_rawReads.txt", "w")
for num in readsHist_T10:
    T10_rawReads.write(str(num))
    T10_rawReads.write("\n")
T10_rawReads.close()
print "make sure nothing changed: ", len(T10_dict), "\n"


print "T11"
#T11 = open("/Users/cindyyeh/Google_Drive/dunham_lab_stuff/allelelibraries/seq_results/61_T11_S5_R1_001.fastq", "r")
timepoint = "T11"
sampleFileName = readFolder+sample+"_"+timepoint+"_merged.assembled.fastq"
T11 = open(sampleFileName,"r")
T11_dict = {}
bc_counts(T11, T11_dict)  # count # of reads per barcode
T11_readct = {}
for key in T11_dict:
	T11_readct[key] = T11_dict[key]
unique_bcs(T11_dict)  # get an idea of # of unique BCs and total # of reads
T11.close()
T11_reads = sum(T11_dict.values())  # total number of reads
readsHist_T11 = []
T11_dict_filtered = {}
normalize(T11_dict, T11_reads, T0_dict, readsHist_T11)  # all values should be normalized now?
# write values to file
T11_hist_out = open(outputLoc+"T11_hist.txt", "w")
writehist(T11_dict, T11_hist_out) # this is for normalized reads
# writehist(T11_dict_filtered, T11_hist_out) # this is for raw read counts
T11_hist_out.close()
T11_rawReads = open(outputLoc+"T11_rawReads.txt", "w")
for num in readsHist_T11:
    T11_rawReads.write(str(num))
    T11_rawReads.write("\n")
T11_rawReads.close()
print "make sure nothing changed: ", len(T11_dict), "\n"


timePointsAll = {}
timePointsToT10 = {}
timePointsToT9 = {}
timePointsToT8 = {}
timePointsToT7 = {}
timePointsToT6 = {}
timePointsToT5 = {}
timePointsToT4 = {}
timePointsFromT2 = {}
timePointsFromT3 = {}
timePointsFromT4 = {}
timePointsFromT5 = {}
timePointsFromT6 = {}
timePointsFromT7 = {}
timePointsFromT8 = {}

for key in T0_dict:
    # make dictionary of barcodes that are present in all 11 timepoints. dict = {entry1 : [freq @ T1, freq @ T2, etc...], entry 2 : [ freq list ]}
    # for each elif statement, it'll remove the latest timepoint that the BC was not sequenced in
    if key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key], T11_dict[key]]
        timePointsAll[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key]]
        timePointsToT10[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key]]
        timePointsToT9[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key]]
        timePointsToT8[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key]]
        timePointsToT7[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key]]
        timePointsToT6[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key]]
        timePointsToT5[key] = tempList
    elif key in T1_dict and key in T2_dict and key in T3_dict and key in T4_dict:
        tempList = [T1_dict[key], T2_dict[key], T3_dict[key], T4_dict[key]]
        timePointsToT4[key] = tempList
    elif key in T2_dict and key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,T2_dict[key], T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT2[key] = tempList
    elif key in T3_dict and key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,0, T3_dict[key], T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT3[key] = tempList
    elif key in T4_dict and key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,0,0, T4_dict[key], T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT4[key] = tempList
    elif key in T5_dict and key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,0,0,0, T5_dict[key], T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT5[key] = tempList
    elif key in T6_dict and key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,0,0,0,0, T6_dict[key], T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT6[key] = tempList
    elif key in T7_dict and key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,0,0,0,0,0, T7_dict[key], T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT7[key] = tempList
    elif key in T8_dict and key in T9_dict and key in T10_dict and key in T11_dict:
        tempList = [0,0,0,0,0,0,0, T8_dict[key], T9_dict[key], T10_dict[key],T11_dict[key]]
        timePointsFromT8[key] = tempList

print "Num of barcodes in entire competition: ", len(timePointsAll)



raw_read_file = open(sample+"_all_read_cts.csv","w+")
#raw_read_file.write("barcode\tT0\tT1\tT2\tT3\tT4\tT5\tT6\tT7\tT9\tT10\tT11")

for key in T0_readct.keys():
	raw_read_file.write("\n"+key)
	if key in T1_readct.keys():
		raw_read_file.write(","+str(T1_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T2_readct:
		raw_read_file.write(","+str(T2_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T3_readct:
		raw_read_file.write(","+str(T3_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T4_readct:
		raw_read_file.write(","+str(T4_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T5_readct:
		raw_read_file.write(","+str(T5_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T6_readct:
		raw_read_file.write(","+str(T6_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T7_readct:
		raw_read_file.write(","+str(T7_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T8_readct:
		raw_read_file.write(","+str(T8_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T9_readct:
		raw_read_file.write(","+str(T9_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T10_readct:
		raw_read_file.write(","+str(T10_readct[key]))
	else:
		raw_read_file.write(","+"0")
	if key in T11_readct:
		raw_read_file.write(","+str(T11_readct[key]))
	else:
		raw_read_file.write(","+"0")
raw_read_file.close()



# this makes a file of all the barcode freqs
fileTimes = open(cwd+"timePoints.txt", "w")
writePoints(timePointsAll, fileTimes)
writePoints(timePointsToT10, fileTimes)
writePoints(timePointsToT9, fileTimes)
writePoints(timePointsToT8, fileTimes)
writePoints(timePointsToT7, fileTimes)
writePoints(timePointsToT6, fileTimes)
writePoints(timePointsToT5, fileTimes)
writePoints(timePointsToT4, fileTimes)
writePoints(timePointsFromT2, fileTimes)
writePoints(timePointsFromT3, fileTimes)
writePoints(timePointsFromT4, fileTimes)
writePoints(timePointsFromT5, fileTimes)
writePoints(timePointsFromT6, fileTimes)
writePoints(timePointsFromT7, fileTimes)
writePoints(timePointsFromT8, fileTimes)
fileTimes.close()
