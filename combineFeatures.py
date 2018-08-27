# Date Created: 22 August 2018
# Date Last Modified: 22 August 2018
# Execution: python combineFeatures.py onlyExpressionMLTableFilename noExpressionMLTableFilename
# argv1: filename for file containing ML table that contains only gene expression features
# argv2: filename for file containing ML table that lacks gene expression features
# Description: TODO
# Run Time: TODO

#!/usr/bin/env python
import sys

# read in command line arguments
onlyExpressionMLTableFilename = sys.argv[1]
noExpressionMLTableFilename = sys.argv[2]

# read in expression ML table
onlyExpressionMLTableFile = open(onlyExpressionMLTableFilename, 'r')

# store header labels
expressionHeaderLine = onlyExpressionMLTableFile.readline()
expressionHeaderLine = expressionHeaderLine.rstrip('\r\n')
expressionHeaders = expressionHeaderLine.split('\t')

numExpressionHeaders = len(expressionHeaders) # for checking
numTissues = len(expressionHeaders) - 2 # doesn't count first two columns, which are snpGroup and snpType

# store each vector in dictionary (key = snp group, value = vector)
onlyExpressionDict = {}

for line in onlyExpressionMLTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	vector = []

	snpGroup = columns[0]
	snpType = columns[1]

	for i in range(numExpressionHeaders):
		vector.append(columns[i])

	onlyExpressionDict[snpGroup] = vector
onlyExpressionMLTableFile.close()

# read in ML table that lacks expression features
noExpressionMLTableFile = open(noExpressionMLTableFilename, 'r')

# store header labels
noExpressionHeaderLine = noExpressionMLTableFile.readline()
noExpressionHeaderLine = noExpressionHeaderLine.rstrip('\r\n')
noExpressionHeaders = noExpressionHeaderLine.split('\t')

numHeaders_noExpression = len(noExpressionHeaders)

# store each vetor in dictionary
noExpressionDict = {}

for line in noExpressionMLTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	vector = []

	snpGroup = columns[0]
	snpType = columns[1]

	for i in range(numHeaders_noExpression):
		vector.append(columns[i])

	noExpressionDict[snpGroup] = vector
noExpressionMLTableFile.close()

# create dictionary to store vector with both features combined
combinedFeaturesDict = {}

for snp in noExpressionDict:
	noExpressionVector = noExpressionDict[snp]
	onlyExpressionVector = onlyExpressionDict[snp]

	# check to make sure snpType is the same --> snp type is in second column (snp group is in first column - index 0)
	if noExpressionDict[snp][1] != onlyExpressionDict[snp][1]:
		print "ERROR: snp types aren't the same!"

	# combine these vectors 
	combinedVector = noExpressionVector

	# remove first two columns of onlyExpressionVector when appending
	for i in range(2, len(onlyExpressionVector)):
		combinedVector.append(onlyExpressionVector[i])

	combinedFeaturesDict[snp] = combinedVector

# create output file
outFilename = noExpressionMLTableFilename.rstrip("_ML_table.txt") + "_combined_ML_table.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# create new header line by combining previous two header lines
newHeaderLine = noExpressionHeaderLine + tab
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