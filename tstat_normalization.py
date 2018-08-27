# Date Created: 3 July 2018
# Date Last Modified: 16 July 2018
# Execution: python tstat_normalization.py tstatFilename
# argv1: filename for file containing tissue expression t-statistics
# Description: This program rank-normalizes the t statistics in the "GTEx.tstat.tsv" file, which is formatted 
# 	such that each row is a gene, and each column is the t statistic for the expression of that gene in the tissue 
# 	corresponding to that column.
# Run Time: 318 sec (~ 5 min)

#!/usr/bin/env python
import sys

# Input: a list
# Output: a list of items that are present in the list more than once
# Description: finds items that are present more than once in the inputted list 
def findDuplicates(anyList):
	seen = set()
	duplicates = []

	for item in anyList:
		if item in seen and item not in duplicates:
			duplicates.append(item)
		else:
			seen.add(item)

	return duplicates

# Input: a list, an item to be searched within the list (optional)
# Output: If item is not specified, returns dictionary containing each item and its respective number of occurrences in the list.
# If item is specified, returns an integer for the number of occurrences of that item in the list. 
# Description: Determines the number of occurences for items in the list. If a particular item is specified, returns the number
# of occurrences of that item in the list. Otherwise, returns a dictionary containing the number of occurences of every item
# in the list. 
def numOccurrences(anyList, item = None):
	numOccurrences = {}

	for i in range(len(anyList)):
		if anyList[i] not in numOccurrences:
			numOccurrences[anyList[i]] = 1
		else: 
			# item is a duplicate
			numOccurrences[anyList[i]] += 1

	if item == None:
		# no item was specified to be parsed in the list
		return numOccurrences
	else:
		# parse list for the number of occurences of the specified item
		if item not in numOccurrences:
			return 0
	
	return numOccurrences[item]

# Input: a list
# Output: a rank-ordered version of the original list
# Description: Determines the ranks of every item in a list, then returns a rank-ordered list, replacing each item in the
# original list with its rank.
def rankList(anyList):
	# sort the inputted list
	sortedList = sorted(anyList)

	# create dictionary to store each item in the list and its corresponding rank
	ranksDict = {}

	# determine whether there are any duplicated items
	duplicates = findDuplicates(anyList)
	if len(duplicates) == 0:
		# list has no duplicated items
		# ranks are simply the indexes of the items in the sorted list
		for item in anyList:
			ranksDict[item] = sortedList.index(item)
	else:
		# list has duplicated items
		print "List has duplicated item."
		for item in anyList:
			# the ranks of unique items are their indexes in the sorted list
			if numOccurrences(anyList, item) == 1:
				ranksDict[item] = sortedList.index(item)
			else:
				# the rank of non-unique items is the average of their indexes in the sorted list
				numItem = numOccurrences(anyList, item)

				firstIndex = float(sortedList.index(item))

				# determine the average of the indexes in the sorted list
				sumRanksSoFar = firstIndex
				for i in range(1, numItem):
					sumRanksSoFar += firstIndex + i

				averageRank = sumRanksSoFar / numItem

				ranksDict[item] = averageRank

	# create the list replacing the items in the original list with their ranks
	rankedList = []
	for item in anyList:
		rank = ranksDict[item]
		rankedList.append(rank)

	return rankedList

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

# initialize matrix comprised of numTissues lists
matrix = [[] for tissue in range(numTissues)]

# create list for ENSGIDs
ids = []

for line in tstatFile:
	line = line.rstrip('\r\n')
	tissues = line.split('\t')

	# add ENSGID, found in first column, to the list
	ids.append(tissues[0])

	# for every subsequent column (representing another tissue), tissues[i + 1] = tissue expression for that gene
	for i in range(numTissues):
		matrix[i].append(float(tissues[i + 1]))

tstatFile.close()
numGenes = len(matrix[0])

# initialize a new matrix to store the ranks
rankedMatrix = [[0 for gene in range(numGenes)] for tissue in range(numTissues)]

# rank tissue expressions for each tissue
for i in range(numTissues):
	rankedMatrix[i] = rankList(matrix[i])

# create new file containing rank-normalized tstats
outFilename = "normalizedGTEx.tstat.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# write header line onto new file
headerLineOutput = headerLine + newline
outFile.write(headerLineOutput)

for i in range(numGenes):
	output = ids[i] + tab
	for j in range(numTissues):
		if j < (numTissues - 1):
			output += str(rankedMatrix[j][i]) + tab
		else: # j = numTissues - 1
			output += str(rankedMatrix[j][i]) + newline

	outFile.write(output)

outFile.close()