# Date Created: 21 August 2018
# Date Last Modified: 21 August 2018
# Execution: python snpsInModel.py [1] [2] [3] [4]
# argv1: filename of file containing lipid testing snps and their type
# argv2: filename of file containing lipid training snps and their type
# argv3: filename of file containing T2D-like testing snps and their type
# argv4: filename of file containing T2D-like training snps and their type
# Description:
# Run Time:

#!/usr/bin/env python
import sys

# read in command line arguments - four files to concatenate
lipidTestingSnpsFilename = sys.argv[1]
lipidTrainingSnpsFilename = sys.argv[2]
T2DLikeTestingSnpsFilename = sys.argv[3]
T2DLikeTrainingSnpsFilename = sys.argv[4]

# read in lipid testing snps
lipidTestingSnpsFile = open(lipidTestingSnpsFilename, 'r')
lipidTestingSnpsFile.readline() # ignore header line

lipidTestingSnpsDict = {}
for line in lipidTestingSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	info = []
	snpGroup = columns[0]
	snpType = columns[1] # whether snp is index or control
	snpCategory = "lipidTesting"

	info.append(snpGroup)
	info.append(snpType)
	info.append(snpCategory)

	lipidTestingSnpsDict[snpGroup] = info

numLipidTestingSnps = len(lipidTestingSnpsDict)
print "Finished reading in the", lipidTestingSnpsFilename, "file, which contained", numLipidTestingSnps, "lipid testing snps."

# read in lipid training snps
lipidTrainingSnpsFile = open(lipidTrainingSnpsFilename, 'r')
lipidTrainingSnpsFile.readline() # ignore header line

lipidTrainingSnpsDict = {}
for line in lipidTrainingSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	info = []
	snpGroup = columns[0]
	snpType = columns[1] # whether snp is index or control
	snpCategory = "lipidTraining"

	info.append(snpGroup)
	info.append(snpType)
	info.append(snpCategory)

	lipidTrainingSnpsDict[snpGroup] = info

numLipidTrainingSnps = len(lipidTrainingSnpsDict)
print "Finished reading in the", lipidTrainingSnpsFilename, "file, which contained", numLipidTrainingSnps, "lipid training snps."

# read in T2D-like testing snps
T2DLikeTestingSnpsFile = open(T2DLikeTestingSnpsFilename, 'r')
T2DLikeTestingSnpsFile.readline() # ignore header line

T2DLikeTestingSnpsDict = {}
for line in T2DLikeTestingSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	info = []
	snpGroup = columns[0]
	snpType = columns[1] # whether snp is index or control
	snpCategory = "T2DLikeTesting"

	info.append(snpGroup)
	info.append(snpType)
	info.append(snpCategory)

	T2DLikeTestingSnpsDict[snpGroup] = info

numT2DLikeTestingSnps = len(T2DLikeTestingSnpsDict)
print "Finished reading in the", T2DLikeTestingSnpsFilename, "file, which contained", numT2DLikeTestingSnps, "T2D-like testing snps."

# read in T2D-like training snps
T2DLikeTrainingSnpsFile = open(T2DLikeTrainingSnpsFilename, 'r')
T2DLikeTrainingSnpsFile.readline() # ignore header line

T2DLikeTrainingSnpsDict = {}
for line in T2DLikeTrainingSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	info = []
	snpGroup = columns[0]
	snpType = columns[1] # whether snp is index or control
	snpCategory = "T2DLikeTraining"

	info.append(snpGroup)
	info.append(snpType)
	info.append(snpCategory)

	T2DLikeTrainingSnpsDict[snpGroup] = info

numT2DLikeTrainingSnps = len(T2DLikeTrainingSnpsDict)
print "Finished reading in the", T2DLikeTrainingSnpsFilename, "file, which contained", numT2DLikeTrainingSnps, "T2D-like training snps."

totalSnps = numLipidTestingSnps + numLipidTrainingSnps + numT2DLikeTestingSnps + numT2DLikeTrainingSnps
print "There are", totalSnps, "being written to the output file."

# create output file
outFilename = "snpTypes.txt"
outFile = open(outFilename, 'w')

tab = '\t'
newline = '\n'

output = "snp" + tab + "type" + tab + "category" + newline # write out the header line

# add lipid testing snps to output file
for snp in lipidTestingSnpsDict:
	for i in range(3): # the info vector contains 3 items for each snp
		if i < 2:
			output += lipidTestingSnpsDict[snp][i] + tab
		else: # add newline at end so each snp is on its own line
			output += lipidTestingSnpsDict[snp][i] + newline
		
# add lipid training snps to output file
for snp in lipidTrainingSnpsDict:
	for i in range(3):
		if i < 2:
			output += lipidTrainingSnpsDict[snp][i] + tab
		else:
			output += lipidTrainingSnpsDict[snp][i] + newline

# add T2D-like testing snps to output file
for snp in T2DLikeTestingSnpsDict:
	for i in range(3):
		if i < 2:
			output += T2DLikeTestingSnpsDict[snp][i] + tab
		else:
			output += T2DLikeTestingSnpsDict[snp][i] + newline

# add T2D-like training snps to output file
for snp in T2DLikeTrainingSnpsDict:
	for i in range(3):
		if i < 2:
			output += T2DLikeTrainingSnpsDict[snp][i] + tab
		else:
			output += T2DLikeTrainingSnpsDict[snp][i] + newline

outFile.write(output)
outFile.close()