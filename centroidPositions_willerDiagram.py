# Date Created: 14 September 2018
# Date Last Modified: 14 September 2018
# Execution: python centroidPositions.py groupedSnpsFilename
# Description: TODO

#!/usr/bin/env python
import sys

# read in command-line arguments
groupedSnpsFilename = sys.argv[1]

# create dictionary that will store groups as keys and the snps in each group as values
snpGroupsDict = {}

# parse groupedSnpsFile
groupedSnpsFile = open(groupedSnpsFilename, 'r')
numSnps = 0
for line in groupedSnpsFile:
	numSnps += 1
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	groupNumber = columns[0]
	snp = columns[1]

	# create list of snps in each group
	snpsInEachGroup = []

	if groupNumber in snpGroupsDict: # group is already stored in the dictionary
		# add the snp to the list of snps contained in that group
		snpGroupsDict[groupNumber].append(snp)
	else: # group isn't already in the dictionary
		# add snp to the list
		snpsInEachGroup.append(snp)
		snpGroupsDict[groupNumber] = snpsInEachGroup
groupedSnpsFile.close()

numSnpGroups = len(snpGroupsDict)
print "There are", numSnpGroups, "groups."

# calculate average number of snps in each group
numSnpsInGroup = []
sumSnps = 0
for i in range(numSnpGroups):
	numSnpsInGroup.append(len(snpGroupsDict.values()[i]))
	sumSnps += len(snpGroupsDict.values()[i])

maxSnps = max(numSnpsInGroup)
minSnps = min(numSnpsInGroup)

averageNumSnpsInGroup = sumSnps / numSnpGroups
print "The number of snps in each group ranges from", minSnps, "to", maxSnps, "with an average of", averageNumSnpsInGroup, "snps." 

# create dictionary that will store the group number as the key and the centroid snp as its value
centroidSnps = {}

for group in snpGroupsDict:
	groupSnpLocations = []

	for i in range(len(snpGroupsDict[group])):
		groupedSnp = snpGroupsDict[group][i]

		# process the snp
		groupedSnp = groupedSnp.split(':')
		chromosome = groupedSnp[0]
		snpLocation = int(groupedSnp[1])

		groupSnpLocations.append(snpLocation)

	minSnpLocation = min(groupSnpLocations)
	maxSnpLocation = max(groupSnpLocations)
	centroidSnpLocation = (minSnpLocation + maxSnpLocation) / 2.0

	centroidSnps[group] = int(centroidSnpLocation)

tab = "\t"
newline = "\n"

output = ""

outFilename = "centroidSnps_" + groupedSnpsFilename
outFile = open(outFilename, 'w')

for group in centroidSnps:
	# obtain chromosome number
	groupChromosome = group.split("_")[0]
	chromosomeNumber = groupChromosome.strip("Chr")

	centroidSnp = "chr" + chromosomeNumber + ":" + str(centroidSnps[group])

	output += group + tab + centroidSnp + newline

outFile.write(output)