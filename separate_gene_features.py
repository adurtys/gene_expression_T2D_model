# Date Created: 24 September 2018
# Date Last Modified: 28 September 2018
# argv[1] = nearestGenes_filename
# argv[2] = normalizedTstat_filename
# argv[3] = groupedSnpTypes_filename
# argv[4] = threshold #

#!/usr/bin/env python
import sys

# read in command-line arguments
nearestGeneMatrixFilename = sys.argv[1]
tstatFilename = sys.argv[2]
groupedSnpTypesFilename = sys.argv[3]
threshold = float(sys.argv[4])

# read in groupedSnpType file
groupedSnpTypesFile = open(groupedSnpTypesFilename, 'r')

# create dictionary containing snp type for each group of snps --> key = snp group, value = vector containing snp group, snp type, and snp source
snpTypeDict = {}

# create dictionary for snps in multiple categories
snpsWithMultipleCategories = {}

# ignore header line in groupedSnpTypesFilename
groupedSnpTypesFile.readline()

numSnps = 0
for line in groupedSnpTypesFile:
	numSnps += 1
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snpGroup = columns[0]
	snpType = columns[1] # "index" or "control"
	snpCategory = columns[2]

	info = []
	info.append(snpGroup)
	info.append(snpType)
	info.append(snpCategory)

	if snpGroup in snpTypeDict:
		# snp has more than one category
		if snpGroup in snpsWithMultipleCategories:
			# snp has more than two categories --> not allowed
			print "ERROR (separate_gene_features.py line 48): snpGroup has more than two categories."
		else:
			snpsWithMultipleCategories[snpGroup] = info
	else:
		snpTypeDict[snpGroup] = info

print "Finished reading in the", groupedSnpTypesFilename, "file, which contained", numSnps, "snps."
groupedSnpTypesFile.close()

# read in normalized t-statistics filename
tstatFile = open(tstatFilename, 'r')

# store column lables in header row of t-statistics file
headerLine = tstatFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

numTissues = len(headers) - 1 # subtract 1 because first column contains gene ID

# make dictionary for ranks of t-statistics --> key = geneId; value = list of tissue-expression ranks
expressionRanks = {}

for line in tstatFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')
	geneId = columns[0]

	tstats = []

	for i in range(numTissues):
		tstats.append(int(columns[i + 1]))

	expressionRanks[geneId] = tstats

tstatFile.close()

genesWithTstats = len(expressionRanks)
print "Finished reading in normalized t-statistics file, which contained tissue expression t-statistics for", genesWithTstats, "genes."

# determine whether expression rank meets threshold for high expression
numTopRankingGenes = threshold * genesWithTstats
critRank = genesWithTstats - numTopRankingGenes

# read in file containing nearest genes
nearestGeneMatrixFile = open(nearestGeneMatrixFilename, 'r')

# make dictionary for nearest genes --> key = snpGroupName; value = list of snpGroupName, centroidSnp, and three nearest genes 
nearestGeneDict = {}
for line in nearestGeneMatrixFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	groupName = columns[0]
	# columns[1] = snp
	# columns[2] = first_nearestGene
	# columns[3] = second_nearestGene
	# columns[4] = third_nearestGene

	nearestGenes = []

	for i in range(5):
		nearestGenes.append(columns[i])

	nearestGeneDict[groupName] = nearestGenes

nearestGeneMatrixFile.close()

# make dictionary that will store gene expression output vectors for the nearest genes --> key = snpGroupName; 
# value = list of snp, snpInfo, 3 output vectors, corresponding to the output vectors for each of the three nearest genes
snpVectorDict = {}

for group in nearestGeneDict:
	# store relevant info for each grouped snp that is in snpTypeDict
	vector_firstNG = []
	vector_secondNG = []
	vector_thirdNG = []
	vector_mergedNG = []
	
	if group in snpTypeDict:
		# centroidSnp = nearestGeneDict[group][1]
		snpType = snpTypeDict[group][1]
		snpCategory = snpTypeDict[group][2]
		first_ng = nearestGeneDict[group][2]
		second_ng = nearestGeneDict[group][3]
		third_ng = nearestGeneDict[group][4]

		# determine output vector for first_nearestGene
		if first_ng not in expressionRanks:
			print "ERROR (separate_gene_features.py line 136): first nearest gene doesn't have t-statistics."
		else:
			for i in range(numTissues):
				if expressionRanks[first_ng][i] >= critRank:
					vector_firstNG.append(1)
				else:
					vector_firstNG.append(0)

		# determine output vector for second_nearestGene
		if second_ng not in expressionRanks:
			print "ERROR (separate_gene_features.py line 146): second nearest gene doesn't have t-statistics."
		else:
			for i in range(numTissues):
				if expressionRanks[second_ng][i] >= critRank:
					vector_secondNG.append(1)
				else:
					vector_secondNG.append(0)

		# determine output vector for third_nearestGene
		if third_ng not in expressionRanks:
			print "ERROR (separate_gene_features.py line 156): third nearest gene doesn't have t-statistics."
		else:
			for i in range(numTissues):
				if expressionRanks[third_ng][i] >= critRank:
					vector_thirdNG.append(1)
				else:
					vector_thirdNG.append(0)

		# determine output vector when merging info for nearest genes
		for i in range(numTissues):
			if (vector_firstNG[i] == 1) or (vector_secondNG[i] == 1) or (vector_thirdNG[i] == 1):
				vector_mergedNG.append(1)
			else:
				vector_mergedNG.append(0)

		# concatenate the four vectors (vector_firstNG, vector_second_NG, vector_thirdNG, vector_mergedNG) into one long vector
		vector = vector_firstNG + vector_secondNG + vector_thirdNG + vector_mergedNG
		vector.insert(0, group)
		vector.insert(1, snpType)

		snpVectorDict[group] = [snpCategory, vector]

##### create the expression ML table ####

# create data structures for the different snp types --> key = snpGroup, value = expression vectors for the three nearest genes and for the combined
lipidTestingSnps = {}
lipidTrainingSnps = {}
T2DLikeTestingSnps = {}
T2DLikeTrainingSnps = {}

for group in snpVectorDict:
	snpCategory = snpVectorDict[group][0]
	snpVector = snpVectorDict[group][1]

	# store snpGroup's info in the appropriate dictionary
	if (snpCategory == "lipidTesting"):
		lipidTestingSnps[group] = snpVector
	elif (snpCategory == "lipidTraining"):
		lipidTrainingSnps[group] = snpVector
	elif (snpCategory == "T2DLikeTesting"):
		T2DLikeTestingSnps[group] = snpVector
	elif (snpCategory == "T2DLikeTraining"):
		T2DLikeTrainingSnps[group] = snpVector
	else:
		print "ERROR: invalid snpCategory!"

	# add snps with multiple categories to both dictionaries
	if group in snpsWithMultipleCategories:
		secondSnpCategory = snpsWithMultipleCategories[snp][2]

		# snp category determines the dictionary in which the snp will be stored
		if (secondSnpCategory == "lipidTesting"):
			lipidTestingSnps[snp] = snpVector
		elif (secondSnpCategory == "lipidTraining"):
			lipidTrainingSnps[snp] = snpVector
		elif (secondSnpCategory == "T2DLikeTesting"):
			T2DLikeTestingSnps[snp] = snpVector
		elif (secondSnpCategory == "T2DLikeTraining"):
			T2DLikeTrainingSnps[snp] = snpVector
		else:
			print "ERROR: invalid snp category!"

totalNumSnps = len(lipidTestingSnps) + len(lipidTrainingSnps) + len(T2DLikeTestingSnps) + len(T2DLikeTrainingSnps)

print "There were", totalNumSnps, "snps in total."
print "There were", len(lipidTestingSnps), "lipid testing snps."
print "There were", len(lipidTrainingSnps), "lipid training snps."
print "There were", len(T2DLikeTestingSnps), "T2D testing snps."
print "There were", len(T2DLikeTrainingSnps), "T2D training snps."
print "There were", len(snpsWithMultipleCategories), "snps in multiple categories."

if numSnps != totalNumSnps:
	print "ERROR: some snps are being lost." # TODO: deal with this!

print "Creating an output file for each of the four types of snp."
# create output files
lipidTestingOutFilename = "lipid_testing_onlyExpression_ML_table.txt"
lipidTrainingOutFilename = "lipid_training_onlyExpression_ML_table.txt"
T2DLikeTestingOutFilename = "T2D_testing_onlyExpression_ML_table.txt"
T2DLikeTrainingOutFilename = "T2D_training_onlyExpression_ML_table.txt"

# open the four output files
lipidTestingOutFile = open(lipidTestingOutFilename, 'w')
lipidTrainingOutFile = open(lipidTrainingOutFilename, 'w')
T2DLikeTestingOutFile = open(T2DLikeTestingOutFilename, 'w')
T2DLikeTrainingOutFile = open(T2DLikeTrainingOutFilename, 'w')

tab = "\t"
newline = "\n"

# create new header line
newHeaderLine = "snp" + tab + "type" + tab

numNewLabels = numTissues * 4
firstNG_indexes = numTissues * 2
secondNG_indexes = numTissues * 3

origLabels = []

for i in range(numNewLabels):
	if i in range(numTissues):
		origLabels.append(headers[i + 1])
	else:
		j = i % numTissues
		if i in range(numTissues, firstNG_indexes):
			newLabel = origLabels[j] + "_firstNG"
			newHeaderLine += newLabel + tab
		elif i in range(firstNG_indexes, secondNG_indexes):
			newLabel = origLabels[j] + "_secondNG"
			newHeaderLine += newLabel + tab
		elif i in range(secondNG_indexes, numNewLabels):
			newLabel = origLabels[j] + "_thirdNG"
			newHeaderLine += newLabel + tab

for i in range(numTissues):
	newLabel = origLabels[i] + "_mergedNG"
	if i < (numTissues - 1):
		newHeaderLine += newLabel + tab
	else: # create new line at end
		newHeaderLine += newLabel + newline

# create lipid testing output file
lipidTestingOutput = newHeaderLine
for group in lipidTestingSnps:
	snpGroup = lipidTestingSnps[group][0]
	snpType = lipidTestingSnps[group][1]

	lipidTestingOutput += snpGroup + tab + snpType + tab

	for i in range(numNewLabels):
		if i < (numNewLabels - 1):
			lipidTestingOutput += str(lipidTestingSnps[group][i + 2]) + tab
		else: # create new line at the end of each vector
			lipidTestingOutput += str(lipidTestingSnps[group][i + 2]) + newline

lipidTestingOutFile.write(lipidTestingOutput)
lipidTestingOutFile.close()

# create lipid training output file
lipidTrainingOutput = newHeaderLine
for group in lipidTrainingSnps:
	snpGroup = lipidTrainingSnps[group][0]
	snpType = lipidTrainingSnps[group][1]

	lipidTrainingOutput += snpGroup + tab + snpType + tab

	for i in range(numNewLabels):
		if i < (numNewLabels - 1):
			lipidTrainingOutput += str(lipidTrainingSnps[group][i + 2]) + tab
		else: # create new line at the end of each vector
			lipidTrainingOutput += str(lipidTrainingSnps[group][i + 2]) + newline

lipidTrainingOutFile.write(lipidTrainingOutput)
lipidTrainingOutFile.close()

# create T2D testing output file
T2DLikeTestingOutput = newHeaderLine
for group in T2DLikeTestingSnps:
	snpGroup = T2DLikeTestingSnps[group][0]
	snpType = T2DLikeTestingSnps[group][1]

	T2DLikeTestingOutput += snpGroup + tab + snpType + tab

	for i in range(numNewLabels):
		if i < (numNewLabels - 1):
			T2DLikeTestingOutput += str(T2DLikeTestingSnps[group][i + 2]) + tab
		else: # create new line at the end of each vector
			T2DLikeTestingOutput += str(T2DLikeTestingSnps[group][i + 2]) + newline

T2DLikeTestingOutFile.write(T2DLikeTestingOutput)
T2DLikeTestingOutFile.close()

# create T2D-like training output file
T2DLikeTrainingOutput = newHeaderLine
for group in T2DLikeTrainingSnps:
	snpGroup = T2DLikeTrainingSnps[group][0]
	snpType = T2DLikeTrainingSnps[group][1]

	T2DLikeTrainingOutput += snpGroup + tab + snpType + tab

	for i in range(numNewLabels):
		if i < (numNewLabels - 1):
			T2DLikeTrainingOutput += str(T2DLikeTrainingSnps[group][i + 2]) + tab
		else: # create new line at the end of each vector
			T2DLikeTrainingOutput += str(T2DLikeTrainingSnps[group][i + 2]) + newline

T2DLikeTrainingOutFile.write(T2DLikeTrainingOutput)
T2DLikeTrainingOutFile.close()

expressionVectorTableFile.close()
print "Finished creating output files."