# Date Created: 25 July 2018
# Date Last Modified: 20 August 2018
# Execution: python expression_lookup.py snpFilename geneAnnotationsFilename tstatFilename numGenes distanceFromSnp
#	threshold outFilename nearestGenesFilename processMissingSnps
# argv1: filename for file containing centroid snps to search (there are two columns in this file - the first is the group number
#	for each centroid snp, the second is the snp itself; each snp is listed on its own line in the format "chr#:location")
# argv2: gene annotations filename
# argv3: filename for file containing rank-orders of the tissue expression t-statistics
# argv4: number of nearest genes to the snp for which to conduct the tissue expression lookup
# argv5: distance from the snp for which to look for nearby genes (in kilobp)
# argv6: the percent of top ranks of t-statistic that should be considered "highly expressed" for each tissue (expressed as a 
#	decimal - ex: 0.10, 0.05)
# argv7: filename for output file containing tissue expression vectors for nearest gene(s) to each snp
# argv8: filename for output file containing nearest gene(s) to each snp
# argv9: flag for how to process snps whose nearest gene has missing tissue expression data, either because there is no
#	gene within the specified distance from the snp or beacuse the nearest gene(s) don't have tissue expression t-statistics
#	if processMissingSnps == "I": find the next nearest gene to the snp that contains tissue expression t-statistics
#	if processMissingSnps == "Z": the tissue expression vector should contain 0 for every tissue
#	if processMissingSnps == "H": do not include snp in final output file of tissue expression vectors
# Description: finds the nearest numGenes to each snp in snpFilename, then assesses whether or not this gene's tissue expression
# 	t-statistics meet the threshold for high expression in each tissue. Outputs a vector of expression in each tissue for every snp,
#	where 1 indicates the nearest gene to the snp is meets the threshold for high expression in a particular tissue. If there is
#	no gene within the specified distance of the snp, or if the nearest gene does not have corresponding tissue-expression
#	t-statistics, these snps are processed as specified in the processMissingSnps flag. Also outputs a file containing each snp
#	and its corresponding nearest numGenes that were used when creating the expression vector.
# Runtime: < 1 second per snp

#!/usr/bin/env python
import sys, bisect, nonOverlappingGenes
"""
Listed below are the functions used in the script.

"""
# Input: a = sorted list (ascending order); x = item to search for in sorted list; n = number of nearest values to x to identify; 
#	side = whether values returned should be to the right or left of x in a
# Output: list of n nearest values to x in a
# Description: function returns a list of n nearest values to x in sorted list a
def nearestNumbers(a, x, n, side):
	nearestIndexes = []

	if x in a:
		nearestIndexes.append(a.index(x))

	if side == "right": # look for values greater than x
		i = bisect.bisect_left(a, x)
		while len(nearestIndexes) < n:
			if i not in nearestIndexes: # x not in a
				nearestIndexes.append(i)
				i += 1
			else: # x in a
				i += 1
				nearestIndexes.append(i)
	elif side == "left": # look for values less than x
		i = bisect.bisect_left(a, x)
		while len(nearestIndexes) < n:
			if i not in nearestIndexes: # x not in a
				nearestIndexes.append(i - 1)
				i -= 1
			else: # x in a
				i -= 1
				nearestIndexes.append(i)
	else:
		print "ERROR: Invalid input value for side. Two options are 'right' and 'left'."
	print nearestIndexes

	# create list of values in sorted array corresponding to each index
	nearestValues = []
	for j in range(len(nearestIndexes)):
		index = nearestIndexes[j]
		if index < len(a):
			nearestValues.append(a[index])

	return nearestValues

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

# MAIN SCRIPT
# error-check for correct number of command-line arguments
print "Arguments passed into expression_lookup3.py:", sys.argv[1:]
if len(sys.argv) != 10:
	print "ERROR (expression_lookup3.py line 95): Incorrect number of command-line arguments!"

# read in command line arguments
snpFilename = sys.argv[1]
geneAnnotationsFilename = sys.argv[2]
tstatFilename = sys.argv[3]
numGenes = int(sys.argv[4])
distanceFromSnp = float(sys.argv[5]) * 1000 # multiply by 1000 to convert kilobp to bp
threshold = float(sys.argv[6])
outFilename = sys.argv[7]
nearestGenesFilename = sys.argv[8]
processMissingSnps = sys.argv[9]

# open the output files
outFile = open(outFilename, 'w')
nearestGenesFile = open(nearestGenesFilename, 'w')

# create data structures to store nearest genes and flagged snps
nearestGenes = {} # key = snp, value = nearest gene(s) (depending on numGenes and processMissingSnps)
noNearbyGene = {} # stores snps that do not have nearby gene(s) within 1mbp on either side (regardless of what processMissingSnps is)
noTissueExpression = {} # stores snps whose nearest gene(s) do not have tissue expression data (regardless of what processMissingSnps is)
equidistantFromGenes = {} # key = snp, value = equidistant genes

# make dictionary for ranks of t-statistics --> key = geneId; value = list of tissue-expression ranks
expressionRanks = {}

# populate this dictionary by parsing file containing rank-ordered t-statistics
tstatFile = open(tstatFilename, 'r')

# store column lables in header row of t-statistics file
headerLine = tstatFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

tab = "\t"
newline = "\n"

# copy all labels from original header line except for the first label (ENSGID)
numLabels = len(headers) - 1
newHeaderLine = ""
for i in range(numLabels):
	if i < (numLabels - 1): # add tab between each header
		newHeaderLine += headers[i + 1] + tab
	else: # i == (numLabels - 1) --> on the last label
		newHeaderLine += headers[i + 1]

numTissues = len(headers) - 1 # subtract 1 because first column contains gene ID

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

# create dictionary containing relevant information pertaining to each gene in the gene annotations file (key = geneId; value = list of info [geneChrom, geneName, geneStart, geneEnd])
geneAnnotations = {}

# populate this dictionary by parsing gene annotations file
geneAnnotationsFile = open(geneAnnotationsFilename, 'r')

for line in geneAnnotationsFile:
	line = line.rstrip('\r\n')
	
	# split line on tab, return list of columns
	columns = line.split('\t')
	geneId = columns[0]
	geneName = columns[1]
	geneChrom = columns[2]
	geneStart = int(columns[3])
	geneEnd = int(columns[4])

	geneInfo = []

	# populate geneInfo with relevant info
	geneInfo.append(geneChrom)
	geneInfo.append(geneName)
	geneInfo.append(geneStart)
	geneInfo.append(geneEnd)

	geneAnnotations[geneId] = geneInfo

geneAnnotationsFile.close()
print "Finished reading in gene annotations file."

# read in genes in the gene annotations file that aren't contained in t-statistics file
genesWithoutTstats = nonOverlappingGenes.getGenesWithoutTstat(tstatFilename, geneAnnotationsFilename)
print "There are", len(genesWithoutTstats), "genes in", geneAnnotationsFilename, "that do not have tissue expression t-statistics."

# create dictionary containing output vectors for each snp --> key = snpName, value = tissue expression vector
snpOutputDict = {}

# parse file containing snps to search
snpFile = open(snpFilename, 'r')

numSnps = 0
for line in snpFile:
	numSnps += 1
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snpGroup = columns[0]
	snp = columns[1]

	print "The centroid snp (to be searched) for", snpGroup, "is:", snp

	# process the snp
	snp = snp.split(':')
	chromosome = snp[0]
	snpLocation = int(snp[1])
	snpName = chromosome + ":" + str(snpLocation)

	# create dictionaries of nearby genes with respect to gene start and end locations
	geneStartLocations = {}
	geneEndLocations = {}

	# parse gene annotations dictionary
	for gene in geneAnnotations:
		# if the gene is on the correct choromsome, add its start and end locations to the dictionaries
		if geneAnnotations[gene][0] == chromosome:
			geneStartLocations[gene] = geneAnnotations[gene][2]
			geneEndLocations[gene] = geneAnnotations[gene][3]

	# sort start and end location dictionaries
	sortedStartLocations = sorted(geneStartLocations.values())
	sortedEndLocations = sorted(geneEndLocations.values())

	# search for nearest genes with respect to start location
	print "Searching for nearest genes with respect to gene start location."
	rightStartLocations = nearestNumbers(sortedStartLocations, snpLocation, numGenes, side="right")
	leftStartLocations = nearestNumbers(sortedStartLocations, snpLocation, numGenes, side="left")
	
	# create list of nearby start locations on either side of the gene
	nearbyStartLocations = leftStartLocations + rightStartLocations

	# create dictionary of nearby genes with respect to start location (key = geneID, value = position)
	genesByStartLocation = {}

	for position in nearbyStartLocations:
		# obtain gene ID from the position value
		geneId = geneStartLocations.keys()[geneStartLocations.values().index(position)]

		if geneId in genesByStartLocation:
			# error check --> overwriting values incorrectly if old position is different from new position for the same gene
			if position != genesByStartLocation[geneId]:
				print "ERROR (expression_lookup3.py line 236):", geneId, "already in genesByStartLocation, with the old position,", genesByStartLocation[geneId], "- we are rewriting this position with the new position,", position
		
		genesByStartLocation[geneId] = position

	print "Nearby genes with respect to gene start locations:", genesByStartLocation

	print "Searching for nearest genes with respect to gene end location."
	rightEndLocations = nearestNumbers(sortedEndLocations, snpLocation, numGenes, side="right")
	leftEndLocations = nearestNumbers(sortedEndLocations, snpLocation, numGenes, side="left")

	# create list of nearby end locations on either side of the gene
	nearbyEndLocations = leftEndLocations + rightEndLocations

	# create list of nearby genes with respect to end location
	genesByEndLocation = {}

	for position in nearbyEndLocations:
		# obtain gene ID from the position value
		geneId = geneEndLocations.keys()[geneEndLocations.values().index(position)]

		if geneId in genesByEndLocation:
			# error check--> in this case, we are overwriting values and something is wrong
			if position != genesByEndLocation[geneId]:
				print "ERROR (expression_lookup3.py line 258):", geneId, "already in genesByEndLocation, with the old position,", genesByEndLocation[geneId], "- we are rewriting this position with the new position,", position

		genesByEndLocation[geneId] = position

	print "Nearby genes with respect to gene end locations:", genesByEndLocation

	# determine closest genes to the snp
	distanceFromSnpDict = {} # key = geneId, value = distance from snp
	for gene in genesByStartLocation:
		for startLocation in nearbyStartLocations:
			if genesByStartLocation[gene] == startLocation: # add gene to dictionary containing distances from the snp
				distanceFromSnpDict[gene] = abs(startLocation - snpLocation)
		print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."
	for gene in genesByEndLocation:
		if gene not in distanceFromSnpDict:
			for endLocation in nearbyEndLocations:
				if genesByEndLocation[gene] == endLocation: # add gene to the dictionary
					distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
		else: # gene is already in distanceFromSnpDict
			# only modify distance if new distance is less than previous distance
			for endLocation in nearbyEndLocations:
				if (genesByEndLocation[gene] == endLocation) and (abs(endLocation - snpLocation) < distanceFromSnpDict[gene]):
					distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
		print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."

	print "Genes closest to the snp:", distanceFromSnpDict
	# sort distanceFromSnpDict by distances
	closestDistances = sorted(distanceFromSnpDict.values())

	# check if there are equidistant genes
	duplicates = findDuplicates(closestDistances)

	# obtain geneIds in the same order as would correspond to the distances after having been sorted in ascending order
	genesToCheck = []
	for i in range(len(closestDistances)):
		sortedGeneId = distanceFromSnpDict.keys()[distanceFromSnpDict.values().index(closestDistances[i])]
		genesToCheck.append(sortedGeneId)

	# determine whether each gene in distanceFromSnpDict should be analyzed for tissue expression
	genesToAnalyze = []

	# create dictionary containing output vectors for each gene being analyzed for a given snp (key = geneId, value = list of outputs (0 if low expression, 1 if high expression))
	geneVectorDict = {}
	geneVector = []

	# default assumption is that the snp is not equidistant from genes, so the expected number of genes to analyze will be the number of genes specified by the user
	expectedNumGenesToAnalyze = numGenes

	if len(duplicates) == 0: # snp is not equidistant from genes
		index = 0
		while len(genesToAnalyze) < expectedNumGenesToAnalyze and (index < len(distanceFromSnpDict)):			
			currentGeneId = genesToCheck[index]
			distanceToCheck = distanceFromSnpDict[currentGeneId]

			# check whether distance of currentGeneId is within specified distance from snp
			if distanceToCheck < distanceFromSnp:
				if currentGeneId not in genesWithoutTstats:
					# gene has tissue expression t-statistics
					genesToAnalyze.append(currentGeneId)
					index += 1
				else: # gene does not have tissue expression t-statistics
					print "(expression_lookup3 line 319): One of the", numGenes, "nearest genes to the snp,", currentGeneId, "doesn't have tissue-expression t-statistic."
					
					# add snp to noTissueExpression (value = currentGeneId)
					noTissueExpression[snpName] = currentGeneId

					# proceed depending on method specified by processMissingSnps
					if processMissingSnps == "I":
						print "Finding next nearest gene within distance and containing tissue-expression t-statistics."
						index += 1
					elif processMissingSnps == "Z":
						print "Tissue expression vector for this gene will be 0 for every tissue."
						
						# add gene to genesToAnalyze
						genesToAnalyze.append(currentGeneId)
						index += 1
					elif processMissingSnps == "H":
						print "Not including this gene in the snp's in final tissue expression output."
						index += 1
			else: # currentGeneId is further from snp than specified by distanceFromSnp
				print "(expression_lookup3 line 338) One of the", numGenes, "nearest genes to the snp,", currentGeneId, "isn't within the distance specified from the snp."
					
				# add snp to noNearbyGene (value = currentGeneId)
				noNearbyGene[snpName] = currentGeneId

				# proceed depending on method specified by processMissingSnps
				if processMissingSnps == "I":
					print "Continuing to assess this gene anyway, because it is one of the nearest", numGenes, "genes to the snp."

					if currentGeneId not in genesWithoutTstats:
						# gene has tissue expression t-statistics
						genesToAnalyze.append(currentGeneId)
						index += 1

					else: # gene does not have tissue expression t-statistics
						print "(expression_lookup3 line 353) One of the", numGenes, "nearest genes to the snp,", currentGeneId, "doesn't have a tissue-expression t-statistic."
						
						# add snp to noTissueExpression
						noTissueExpression[snpName] = currentGeneId

						print "Finding the next nearest gene containing tissue-expression t-statistics."
						index += 1
				elif processMissingSnps == "Z":
					print "Tissue expression for this gene will be 0 for every tissue."

					# add gene to genesToAnalyze
					genesToAnalyze.append(currentGeneId)
					index += 1
				elif processMissingSnps == "H":
					print "Not including this gene in the snp's final tissue expression output."
					index += 1
			
	else: # snp is equidistant from genes --> assumes only equidistant from two genes and from one pair of genes
		print "snp is equidistant from genes."
		
		# add snp to equidistantFromGenes
		equidistantFromGenes[snpName] = duplicates

		# assuming gene is only equidistant from two genes (and from one pair of genes) --> store one more gene than what was specified
		expectedNumGenesToAnalyze = numGenes + 1

		index = 0
		while len(genesToAnalyze) < expectedNumGenesToAnalyze and (index < len(distanceFromSnpDict)):
			currentGeneId = genesToCheck[index]
			distanceToCheck = distanceFromSnpDict[currentGeneId]

			# check whether distance of currentGeneId is within specified distance from snp
			if distanceToCheck < distanceFromSnp:
				if currentGeneId not in genesWithoutTstats:
					# gene has tissue expression t-statistics

					if currentGeneId in duplicates: # snp is equidistant from this gene and the next gene in the list
						nextGeneId = genesToCheck[index + 1]

						print snpName, "is equidistant from", currentGeneId, "and", nextGeneId
						print "Will proceed by analyzing both genes."
							
						genesToAnalyze.append(currentGeneId)
						genesToAnalyze.append(nextGeneId)
						index += 2 # increment index by two because assessed and added both current and subsequent genes

					else: # snp isn't equidistnat from this gene and another --> add only this gene
						genesToAnalyze.append(currentGeneId)
						index += 1

				else: # gene does not have tissue expression t-statistics
					print "(expression_lookup3 line 404): One of the", numGenes, "nearest genes to the snp,", currentGeneId, "doesn't have tissue-expression t-statistic."
						
					# add snp to noTissueExpression (value = currentGeneId)
					noTissueExpression[snpName] = currentGeneId

					# proceed depending on method specified by processMissingSnps
					if processMissingSnps == "I":
						print "Finding next nearest gene within distance and containing tissue-expression t-statistics."
						index += 1
					elif processMissingSnps == "Z":
						print "Tissue expression vector for this gene will be 0 for every tissue."
							
						# add gene to genesToAnalyze
						genesToAnalyze.append(currentGeneId)
						index += 1
					elif processMissingSnps == "H":
						print "Not including this gene in the snp's in final tissue expression output."
						index += 1
			else: # currentGeneId is further from snp than specified by distanceFromSnp
				print "(expression_lookup3 line 423) One of the", numGenes, "nearest genes to the snp,", currentGeneId, "isn't within the distance specified from the snp."
					
				# add snp to noNearbyGene (value = currentGeneId)
				noNearbyGene[snpName] = currentGeneId

				# proceed depending on method specified by processMissingSnps
				if processMissingSnps == "I":
					print "Continuing to assess this gene anyway, because it is one of the nearest", numGenes, "genes to the snp."

					if currentGeneId not in genesWithoutTsats:
						# gene has tissue expression t-statistics
						if currentGeneId in duplicates: # snp is equidistant from this gene and the next gene in the list
							nextGeneId = genesToCheck[index + 1]

							print snpName, "is equidistant from", currentGeneId, "and", nextGeneId
							print "Will proceed by analyzing both genes."
							
							genesToAnalyze.append(currentGeneId)
							genesToAnalyze.append(nextGeneId)
							index += 2 # increment index by two because assessed and added both current and subsequent genes

						else: # snp isn't equidistant from this gene and another --> add only this gene
							genesToAnalyze.append(currentGeneId)
							index += 1

					else: # gene does not have tissue expression t-statistics
						print "(expression_lookup3 line 449) One of the", numGenes, "nearest genes to the snp,", currentGeneId, "doesn't have a tissue-expression t-statistic."
							
						# add snp to noTissueExpression
						noTissueExpression[snpName] = currentGeneId

						print "Finding the next nearest gene containing tissue-expression t-statistics."
						index += 1
				elif processMissingSnps == "Z":
					print "Tissue expression for this gene will be 0 for every tissue."

					# add gene to genesToAnalyze
					genesToAnalyze.append(currentGeneId)
					index += 1
				elif processMissingSnps == "H":
					print "Not including this gene in the snp's final tissue expression output."
					index += 1

	print "There are", len(genesToAnalyze), "genes to analyze for tissue expression:", genesToAnalyze

	if (processMissingSnps == "I") and (len(genesToAnalyze) != expectedNumGenesToAnalyze):
		# couldn't find a gene in distanceFromSnpDict that was within the distance specified and also had tissue expression t-statistics
		print "(expression_lookup3.py line 470): Couldn't find a gene in distanceFromSnpDict that was within the distance specified and also had tissue expression t-statistics. Including a larger distance and more genes in the search."

		index = 0
		newNumGenes = numGenes
		newDistanceFromSnp = distanceFromSnp
		while (index < len(distanceFromSnpDict)) and (len(genesToAnalyze) != expectedNumGenesToAnalyze):
			newDistanceFromSnp *= 2 # increase distance from the snp within which to search for genes
			newNumGenes *= 2 # increase number of genes to find

			# search for nearest genes with respect to start location
			print "Searching for nearest genes with respect to gene start location."
			rightStartLocations = nearestNumbers(sortedStartLocations, snpLocation, newNumGenes, side="right")
			leftStartLocations = nearestNumbers(sortedStartLocations, snpLocation, newNumGenes, side="left")
			
			# create list of nearby start locations on either side of the gene
			nearbyStartLocations = leftStartLocations + rightStartLocations

			# create dictionary of nearby genes with respect to start location (key = geneID, value = position)
			genesByStartLocation = {}

			for position in nearbyStartLocations:
				# obtain gene ID from the position value
				geneId = geneStartLocations.keys()[geneStartLocations.values().index(position)]

				if geneId in genesByStartLocation:
					# error check --> in this case, we are overwriting values and something is wrong
					if position != genesByStartLocation[geneId]:
						print "ERROR (expression_lookup3.py line 497):", geneId, "already in genesByStartLocation, with the old position,", genesByStartLocation[geneId], "- we are rewriting this position with the new position,", position
		
				genesByStartLocation[geneId] = position

			print "Nearby genes with respect to gene start locations:", genesByStartLocation

			print "Searching for nearest genes with respect to gene end location."
			rightEndLocations = nearestNumbers(sortedEndLocations, snpLocation, newNumGenes, side="right")
			leftEndLocations = nearestNumbers(sortedEndLocations, snpLocation, newNumGenes, side="left")

			# create list of nearby end locations on either side of the gene
			nearbyEndLocations = leftEndLocations + rightEndLocations

			# create list of nearby genes with respect to end location
			genesByEndLocation = {}

			for position in nearbyEndLocations:
				# obtain gene ID from the position value
				geneId = geneEndLocations.keys()[geneEndLocations.values().index(position)]

				if geneId in genesByEndLocation:
					# error check--> in this case, we are overwriting values and something is wrong
					if position != genesByEndLocation[geneId]:
						print "ERROR (expression_lookup3.py line 520):", geneId, "already in genesByEndLocation, with the old position,", genesByEndLocation[geneId], "- we are rewriting this position with the new position,", position
		
				genesByEndLocation[geneId] = position

			print "Nearby genes with respect to gene end locations:", genesByEndLocation

			# determine closest genes to the snp
			distanceFromSnpDict = {} # key = geneId, value = distance from snp
			for gene in genesByStartLocation:
				for startLocation in nearbyStartLocations:
					if genesByStartLocation[gene] == startLocation: # add gene to dictionary containing distances from the snp
						distanceFromSnpDict[gene] = abs(startLocation - snpLocation)
				print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."
			for gene in genesByEndLocation:
				if gene not in distanceFromSnpDict:
					for endLocation in nearbyEndLocations:
						if genesByEndLocation[gene] == endLocation: # add gene to the dictionary
							distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
				else: # gene is already in distanceFromSnpDict
					# only modify distance if new distance is less than previous distance
					for endLocation in nearbyEndLocations:
						if (genesByEndLocation[gene] == endLocation) and (abs(endLocation - snpLocation) < distanceFromSnpDict[gene]):
							distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
				print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."

			print "Genes closest to the snp:", distanceFromSnpDict
			# sort distanceFromSnpDict by distances
			closestDistances = sorted(distanceFromSnpDict.values())

			# check if there are equidistant genes
			duplicates = findDuplicates(closestDistances)

			# obtain geneIds in the same order as would correspond to the distances after having been sorted in ascending order
			genesToCheck = []
			for i in range(len(closestDistances)):
				sortedGeneId = distanceFromSnpDict.keys()[distanceFromSnpDict.values().index(closestDistances[i])]
				genesToCheck.append(sortedGeneId)

			# determine whether each gene in distanceFromSnpDict should be analyzed for tissue expression
			genesToAnalyze = []

			# create dictionary containing output vectors for each gene being analyzed for a given snp (key = geneId, value = list of outputs (0 if low expression, 1 if high expression))
			geneVectorDict = {}
			geneVector = []

			# NOTE: didn't modify expectedNumGenesToAnalyze because still only want however many genes was initially inputted by the user

			if len(duplicates) == 0: # snp is not equidistant from genes
				currentGeneId = genesToCheck[index]
				distanceToCheck = distanceFromSnpDict[currentGeneId]

				# check whether distance of currentGeneId is within specified distance from snp
				if distanceToCheck < newDistanceFromSnp:
					if currentGeneId not in genesWithoutTstats:
						# gene has tissue expression t-statistics
						genesToAnalyze.append(currentGeneId)
						index += 1
					else: # gene does not have tissue expression t-statistics
						print "(expression_lookup3 line 578): One of the", newNumGenes, "nearest genes to the snp,", currentGeneId, "doesn't have tissue-expression t-statistic."
							
						# add snp to noTissueExpression (value = currentGeneId)
						noTissueExpression[snpName] = currentGeneId

						print "Finding next nearest gene within distance and containing tissue-expression t-statistics."
						index += 1

				else: # currentGeneId is further from snp than specified by distanceFromSnp
					print "(expression_lookup3 line 587) One of the", newNumGenes, "nearest genes to the snp,", currentGeneId, "isn't within the distance specified from the snp."
					
					# add snp to noNearbyGene (value = currentGeneId)
					noNearbyGene[snpName] = currentGeneId

					print "Continuing to assess this gene anyway, because it is one of the nearest", newNumGenes, "genes to the snp."

					if currentGeneId not in genesWithoutTstats:
						# gene has tissue expression t-statistics
						genesToAnalyze.append(currentGeneId)
						index += 1

					else: # gene does not have tissue expression t-statistics
						print "(expression_lookup3 line 600) One of the", newNumGenes, "nearest genes to the snp,", currentGeneId, "doesn't have a tissue-expression t-statistic."
							
						# add snp to noTissueExpression
						noTissueExpression[snpName] = currentGeneId

						print "Finding the next nearest gene containing tissue-expression t-statistics."
						index += 1
			else: # snp is equidistant from genes --> assumes only equidistant from two genes and from one pair of genes
				print "snp is equidistant from genes."
				
				# add snp to equidistantFromGenes
				equidistantFromGenes[snpName] = duplicates

				# assuming gene is only equidistant from two genes (and from one pair of genes) --> store one more gene than what was specified
				expectedNumGenesToAnalyze = numGenes + 1


				currentGeneId = genesToCheck[index]
				distanceToCheck = distanceFromSnpDict[currentGeneId]

				# check whether distance of currentGeneId is within specified distance from snp
				if distanceToCheck < newDistanceFromSnp:
					if currentGeneId not in genesWithoutTstats:
						# gene has tissue expression t-statistics

						if currentGeneId in duplicates: # snp is equidistant from this gene and the next gene in the list
							nextGeneId = genesToCheck[index + 1]

							print snpName, "is equidistant from", currentGeneId, "and", nextGeneId
							print "Will proceed by analyzing both genes."
									
							genesToAnalyze.append(currentGeneId)
							genesToAnalyze.append(nextGeneId)
							index += 2 # increment index by two because assessed and added both current and subsequent genes

						else: # snp isn't equidistant from this gene and another --> add only this gene
							genesToAnalyze.append(currentGeneId)
							index += 1

					else: # gene does not have tissue expression t-statistics
						print "(expression_lookup3 line 640): One of the", newNumGenes, "nearest genes to the snp,", currentGeneId, "doesn't have tissue-expression t-statistic."
								
						# add snp to noTissueExpression (value = currentGeneId)
						noTissueExpression[snpName] = currentGeneId

						print "Finding next nearest gene within distance and containing tissue-expression t-statistics."
						index += 1

				else: # currentGeneId is further from snp than specified by distanceFromSnp
					print "(expression_lookup3 line 649) One of the", newNumGenes, "nearest genes to the snp,", currentGeneId, "isn't within the distance specified from the snp."
							
					# add snp to noNearbyGene (value = currentGeneId)
					noNearbyGene[snpName] = currentGeneId

					print "Continuing to assess this gene anyway, because it is one of the nearest", newNumGenes, "genes to the snp."

					if currentGeneId not in genesWithoutTsats:
						# gene has tissue expression t-statistics
						if currentGeneId in duplicates: # snp is equidistant from this gene and the next gene in the list
							nextGeneId = genesToCheck[index + 1]

							print snpName, "is equidistant from", currentGeneId, "and", nextGeneId
							print "Will proceed by analyzing both genes."
								
							genesToAnalyze.append(currentGeneId)
							genesToAnalyze.append(nextGeneId)
							index += 2 # increment index by two because assessed and added both current and subsequent genes

						else: # snp isn't equidistant from this gene and another --> add only this gene
							genesToAnalyze.append(currentGeneId)
							index += 1

					else: # gene does not have tissue expression t-statistics
						print "(expression_lookup3 line 673) One of the", newNumGenes, "nearest genes to the snp,", currentGeneId, "doesn't have a tissue-expression t-statistic."
								
						# add snp to noTissueExpression
						noTissueExpression[snpName] = currentGeneId

						print "Finding the next nearest gene containing tissue-expression t-statistics."
						index += 1

			print "Upon recreating distanceFromSnpDict, there are", len(genesToAnalyze), "genes to analyze for tissue expression:", genesToAnalyze

	# proceed with analysis
	if len(genesToAnalyze) == expectedNumGenesToAnalyze: 
		# this should be the case for processMissingSnps = I or Z (could also be the case for H, but not necessarily)
		for gene in genesToAnalyze:
			if (gene in noNearbyGene.values()) or (gene in noTissueExpression.values()): # would only keep such a gene if processMissingSnps == Z or I (not if it equals H)
				if processMissingSnps == "H":
					print "ERROR (expression_lookup3.py line 689):", gene, "is being analyzed, but either doesn't have t-statistic or isn't within the specified distance of the snp."
					
				else: # processMissingSnps == Z
					# expression vector is 0 for every tissue for the gene
					for i in range(numTissues):
						geneVector.append(0)

					geneVectorDict[gene] = geneVector

			else:
				# determine whether expression vector meets threshold and add 1's and 0's accordingly
				for i in range(numTissues):
					if expressionRanks[gene][i] >= critRank:
						geneVector.append(1)
					else:
						geneVector.append(0)

				geneVectorDict[gene] = geneVector
	else: # len(genesToAnalyze) != expectedNumGenesToAnalyze (this should only be the case for processMissingSnps = H)
		if processMissingSnps != "H":
			print "ERROR (expression_lookup3.py line 709): Not analyzing the same number of genes as was specified by the user."
			
		if len(genesToAnalyze) != 0: # there are genes to analyze for the snp
			for gene in genesToAnalyze:
				# determine whether expression vector meets threshold and add 1's and 0's accordingly
				for i in range(numTissues):
					if expressionRanks[gene][i] >= critRank:
						geneVector.append(1)
					else:
						geneVector.append(0)

				geneVectorDict[gene] = geneVector

	print "Gene Vector Dictionary:", geneVectorDict

	# only store snp in output dictionary if a gene was stored for the snp
	if len(geneVectorDict) != 0:
		print "Analyzing expression for:", geneVectorDict.keys()
		
		# create tissue expression vector for each snp
		snpOutputVector = []
		
		for i in range(numTissues):
			# initialize this vector as 0 for all tissues
			snpOutputVector.append(0)

		# combine all the geneVectorDict values for each snp	
		for gene in geneVectorDict:
			for i in range(numTissues):
				if geneVectorDict[gene][i] == 1:
					# modify the snp output vector if the tissue expression for any gene in the geneVector dictionary has 
					snpOutputVector[i] = 1

		# use snpOutputVector as the key in the output dictionary
		snpOutputDict[snpName] = snpOutputVector

	else: # don't include snp in the dictionary (or the output file) if no gene expression vector was stored for it
		# len(geneVectorDict == 0) (can only occur if processMissingSnps == H and len(genesToAnalyze) = 0)
		if (processMissingSnps != "H") or len(genesToAnalyze != 0):
			print "ERROR (expression_lookup3.py line 748): Error processing expression vectors."

		print "Excluding", snpName, "from the output file because there are no genes with tissue expression associated with it."

	# create nearestGenes dictionary for every snp, even if it doesn't contain a gene to analyze
	nearestGenes[snpName] = genesToAnalyze

	# create output files
	if numSnps == 1:
		# write header line onto the new output file if the snp being analyzed is the first snp in the list
		headerLineOutput = "snp" + tab + newHeaderLine + newline # edited headerLine so that first column is snp (not ENSGID)
		outFile.write(headerLineOutput)

	# create expression lookup file if the snp is to be included in the output file
	if snpName in snpOutputDict:
		output = snpGroup + tab
		for i in range(numTissues):
			if i < (numTissues - 1):
				output += str(snpOutputDict[snpName][i]) + tab
			else: # add newline if this is the last entry in the snp output vector
				output += str(snpOutputDict[snpName][i]) + newline	
		outFile.write(output)

	# create file containing nearest genes to each snp, determined as specified by processMissingSnps
	nearestGeneOutput = snpGroup + tab + snpName + tab
	for i in range(len(genesToAnalyze)):
		if i < (len(genesToAnalyze) - 1):
			nearestGeneOutput += genesToAnalyze[i] + tab
		else: # add newline if this is the last gene in the list of genes associated with the given snp
			nearestGeneOutput += genesToAnalyze[i] + newline
	nearestGenesFile.write(nearestGeneOutput)

outFile.close()
nearestGenesFile.close()
snpFile.close()

print "There were", numSnps, "snps for which to conduct the tissue expression lookup."
print "There were", len(equidistantFromGenes), "instances in which a snp was equidistant from two genes:", equidistantFromGenes.keys()
print "There were", len(noNearbyGene), "snps which did not have a gene within the distance specified:", noNearbyGene.keys()
print "There were", len(noTissueExpression), "snps which did not have tissue-expression t-statistics for their nearest gene(s):", noTissueExpression.keys()
print "In total, there were", len(genesWithoutTstats), "genes in the geneAnnotationsFile that did not have corresponding tissue-expression t-statistics."
print "These cases of distant genes or missing tissue expression t-statistics were dealt with as specified by the processMissingSnps flag:", processMissingSnps