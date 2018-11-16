# Date Created: 17 October 2018
# Date Last Modified: 17 October 2018
# Execution: python tstat_normalization.py tstatFilename
# argv1: filename for file containing tissue expression t-statistics
# Description: This program normalizes the t-statistics in the "GTEx.tstat.tsv" file, which is formatted 
# 	such that each row is a gene, and each column is the t statistic for the expression of that gene in the tissue 
# 	corresponding to that column. The normalization is done by scaling all the t-statistics from 0 to 1.
# Run Time: ~2 sec

#!/usr/bin/env python
import sys

# check to make sure file was run with correct number of arguments
if len(sys.argv) != 2:
	print "ERROR (tstat_normaliation.py line 104): Incorrect number of command-line arguments!"

# read in the tstat file
tstatFilename = sys.argv[1]
tstatFile = open(tstatFilename, 'r')

# store column labels in header row
headerLine = tstatFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

# exclude first column when counting the number of tissues present in the file
numTissues = len(headers) - 1

# create matrix --> a dictionary whose key = ensgID and values = list of tstats
matrix = {}

for line in tstatFile:
	line = line.rstrip('\r\n')
	tissues = line.split('\t')

	geneId = tissues[0]

	tstats = []
	for i in range(numTissues):
		tstats.append(float(tissues[i + 1]))

	matrix[geneId] = tstats

tstatFile.close()

numGenes = len(matrix)

# initialize a new matrix to store the scaled t-stats --> key = geneId, value = scaled t-stats
scaledMatrix = {}

# scale tissue expression t-stats for each tissue
for gene in matrix:
	maxStat = max(matrix[gene])
	minStat = min(matrix[gene])

	scaledTstats = []
	for i in range(numTissues):
		scaled = (matrix[gene][i] - minStat) / float((maxStat - minStat))
		scaledTstats.append(scaled)

	scaledMatrix[gene] = scaledTstats

# create new file containing rank-normalized tstats
outFilename = "scaledGTEx.tstat.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

output = headerLine + newline

for gene in scaledMatrix:
	output += gene + tab
	for i in range(numTissues):
		if i < (numTissues - 1):
			output += str(scaledMatrix[gene][i]) + tab
		else: # i = numTisses - 1
			output += str(scaledMatrix[gene][i]) + newline

outFile.write(output)
outFile.close()