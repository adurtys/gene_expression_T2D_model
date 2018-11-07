# Date Created: 22 August 2018
# Date Last Modified: 5 September 2018
# Execution: python combineFeatures.py lowExpressionMLTableFilename highExpressionMLTableFilename
# argv1: filename for file containing ML table that contains low gene expression features
# argv2: filename for file containing ML table that lacks gene expression features
# Description: combines gene expression vector to vector of all other regulatory features for all snps contained in highExpressionMLTableFilename
# Run Time: short

#!/usr/bin/env python
import sys

# read in command line arguments
lowExpressionMLTableFilename = sys.argv[1]
highExpressionMLTableFilename = sys.argv[2]

# read in expression ML table
lowExpressionMLTableFile = open(lowExpressionMLTableFilename, 'r')

# store header labels
lowExpressionHeaderLine = lowExpressionMLTableFile.readline()
lowExpressionHeaderLine = lowExpressionHeaderLine.rstrip('\r\n')
lowExpressionHeaders = lowExpressionHeaderLine.split('\t')

numExpressionHeaders = len(lowExpressionHeaders) # for checking
numTissues = len(lowExpressionHeaders) - 2 # doesn't count first two columns, which are snpGroup and snpType

# store each vector in dictionary (key = snp group, value = vector)
lowExpressionDict = {}

for line in lowExpressionMLTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	vector = []

	snpGroup = columns[0]
	snpType = columns[1]

	for i in range(numExpressionHeaders):
		vector.append(columns[i])

	lowExpressionDict[snpGroup] = vector
lowExpressionMLTableFile.close()

# read in ML table that lacks expression features
highExpressionMLTableFile = open(highExpressionMLTableFilename, 'r')

# store header labels
highExpressionHeaderLine = highExpressionMLTableFile.readline()
highExpressionHeaderLine = highExpressionHeaderLine.rstrip('\r\n')
highExpressionHeaders = highExpressionHeaderLine.split('\t')

numHeaders_highExpression = len(highExpressionHeaders)

# store each vetor in dictionary
highExpressionDict = {}

for line in highExpressionMLTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	vector = []

	snpGroup = columns[0]
	snpType = columns[1]

	for i in range(numHeaders_highExpression):
		vector.append(columns[i])

	highExpressionDict[snpGroup] = vector
highExpressionMLTableFile.close()

# create dictionary to store vector with both features combined
combinedFeaturesDict = {}

for snp in highExpressionDict:
	if snp in lowExpressionDict:
		highExpressionVector = highExpressionDict[snp]
		lowExpressionVector = lowExpressionDict[snp]

		# snp type is the one specified in highExpressionDict
		# combine these vectors 
		combinedVector = highExpressionVector

		# remove first two columns of lowExpressionVector when appending
		for i in range(2, len(lowExpressionVector)):
			combinedVector.append(lowExpressionVector[i])

		combinedFeaturesDict[snp] = combinedVector

# create output file TODO: EDDIT THIS!!!!!!!!
outFilename = highExpressionMLTableFilename.rstrip("_highExpressionOnly.txt") + "_combined_lowAndHigh.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# adjust labels in both header lines to indicate whether feature is for high or low expression
newHighExpressionHeaders = []
for header in highExpressionHeaders:
	newHeader = "highExp_" + header
	newHighExpressionHeaders.append(newHeader)

newLowExpressionHeaders = []
for header in lowExpressionHeaders:
	newHeader = "lowExp_" + header
	newLowExpressionHeaders.append(newHeader)

# create new header line by combining two new header lines
newHeaderLine = ""
for header in newHighExpressionHeaders:
	newHeaderLine += header + tab

for i in range(2, numExpressionHeaders):
	if i < (numExpressionHeaders - 1):
		newHeaderLine += expressionHeaders[i] + tab
	else: # add new line at end
		newHeaderLine += expressionHeaders[i] + newline

output = newHeaderLine

for snp in combinedFeaturesDict:
	vector = combinedFeaturesDict[snp]

	for i in range(len(vector)):
		if i < (len(vector) - 1):
			output += vector[i] + tab
		else: # add new line at the end
			output += vector[i] + newline

outFile.write(output)
outFile.close()