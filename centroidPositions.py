# Date Created: 24 July 2018
# Date Last Modified: 3 August 2018
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
for line in groupedSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	groupNumber = int(columns[0])
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
print "There are", len(snpGroupsDict), "groups."

# create dictionary that will store the group number as the key and the centroid snp as its value
centroidSnps = {}

chromosomeNumber = -1
for group in snpGroupsDict:
	groupSnpLocations = []

	for i in range(len(snpGroupsDict[group])):
		groupedSnp = snpGroupsDict[group][i]

		# process the snp
		groupedSnp = groupedSnp.split(':')
		chromosome = groupedSnp[0]
		snpLocation = int(groupedSnp[1])

		groupSnpLocations.append(snpLocation)

		if (len(centroidSnps) == 0) and (i == 0):
			# obtain chromosome number
			chromosomeNumber = int(chromosome.strip("chr"))
			print "These grouped snps are for chromosome", chromosomeNumber

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
	groupName = "Chr" + str(chromosomeNumber) + "_Group_" + str(group)
	centroidSnp = "chr" + str(chromosomeNumber) + ":" + str(centroidSnps[group])

	output += groupName + tab + centroidSnp + newline

outFile.write(output)