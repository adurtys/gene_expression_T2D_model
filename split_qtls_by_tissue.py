# Date Created: 4 May 2019
# Date Last Modified: 4 May 2019
# Execution: python split_qtls_by_tissue.py [1]
# argv[1] = isLD.eQTLs.forTable.sortedByTissue
# Description: split up eQTLs_byTissue_file into separate files for each tissue
# Output: creates directory containing (TODO: FIGURE OUT NUMBER) files (one for each tissue) --> each file contains snpgene pairs in that tissue
# Runtime: TODO

#!/usr/bin/env python
import sys, os

# read in significantQTLs file sorted by tissue
eQTLs_byTissue_filename = sys.argv[1]
eQTLs_byTissue_file = open(eQTLs_byTissue_filename, 'r')

eQTLs_byTissue_dict = {} # key = tissue, value = array of all snpgene_pairs in the tissue (format: eQTL_geneRegulated)
numErrorSnpGenePairs = 0
for line in eQTLs_byTissue_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	eQTL = columns[0]
	geneRegulated = columns[1]
	tissue = columns[2]

	snpgene_pair = eQTL + "_" + geneRegulated
	if tissue not in eQTLs_byTissue_dict:
		eQTLs_byTissue_dict[tissue] = [snpgene_pair]
	else:
		if snpgene_pair not in eQTLs_byTissue_dict[tissue]:
			eQTLs_byTissue_dict[tissue].append(snpgene_pair)
		else:
			print "ERROR: snpgene_pair already stored for", tissue
			numErrorSnpGenePairs += 1

eQTLs_byTissue_file.close()

# create a new directory to store files for each gene
outputDirectoryPath = "./eQTLs_byTissue_directory"
os.mkdir(outputDirectoryPath)

# change directory to output directory
os.chdir(outputDirectoryPath)

# create output files for each tissue
tab = "\t"
newline = "\n"

for tissue in eQTLs_byTissue_dict:
	# create output file
	outputFilename = tissue + "_eQTLs.txt"
	outputFile = open(outputFilename, 'w')

	output = ""
	for snpgene_pair in eQTLs_byTissue_dict[tissue]:
		snpgene_pair = snpgene_pair.split("_")

		output += snpgene_pair[0] + tab + snpgene_pair[1] + newline

	outputFile.close()




