# Date Created: 4 May 2019
# Date Last Modified: 4 May 2019
# Execution: python isQTL_byTissue.py [1] [2] [3]
# argv[1]: directory with eQTLs split by tissue (separate file for each tissue) (/project/voight_ML/adurtys/eqtl_feature/new_2019/isQTL_specificTissue_anyGene/eQTLs_byTissue_directory)
# argv[2]: matrix for eQTL snps that are in LD (r2 > 0.7) w/ GWAS snps (isLD.txt)
# argv[3]: grouped snps file (final_grouped_snps_all.txt)
# Description: TODO

#!/usr/bin/env python
import sys, os

# read in command line arguments
eQTL_byTissue_direc = sys.argv[1]
eQTL_inLD_wGWAS_filename = sys.argv[2]
groupedSnpsFilename = sys.argv[3]

eQTL_tissue_filenames = os.listdir(eQTL_byTissue_direc)

# read in tissue information
eQTL_dict = {} # key = tissue, value = array of significant eQTL snps for each tissue
for filename in eQTL_tissue_filenames:
	tissue = filename.rstrip("_Analysis_eQTLs.txt")
	filepath = eQTL_byTissue_direc + "/" + filename

	file = open(filepath, 'r')

	eQTL_snps_inTissue = []

	for line in file:
		line = line.rstrip('\r\n')
		columns = line.split('\t')

		eQTL = columns[0]
		# geneRegulated = columns[1]

		eQTL = eQTL.split("_")
		eQTL = "chr" + eQTL[0] + ":" + eQTL[1]

		eQTL_snps_inTissue.append(eQTL)

	eQTL_dict[tissue] = eQTL_snps_inTissue
	
	file.close()

# read in eQTL snps in LD w/ GWAS snps
eQTL_inLD_wGWAS_file = open(eQTL_inLD_wGWAS_filename, 'r')

# skip headerline
eQTL_inLD_wGWAS_file.readline()

gwas_isQTL_dict = {} # key = GWAS variantID, value = [GWAS rsID, GWAS variant ID, eQTL rsID, eQTL variant ID]
for line in eQTL_inLD_wGWAS_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')
	
	gwas_rsID = columns[1]
	eQTL_variant_ID = columns[2]
	eQTL_rsID = columns[3]	

	# reformat gwas variant id to match the one used for the snps in the model
	gwas_variantID = columns[0]
	gwas_variantID = gwas_variantID.split("_")
	gwas_variantID = "chr" + gwas_variantID[0] + ":" + gwas_variantID[1]

	gwas_isQTL_dict[gwas_variantID] = [gwas_variantID, gwas_rsID, eQTL_rsID, eQTL_variant_ID]

eQTL_inLD_wGWAS_file.close()

print "There are", len(gwas_isQTL_dict), "GWAS snps in the dictionary."

# read in command-line arguments
groupedSnpsFilename = sys.argv[3]

# create dictionary that will store groups as keys and the snps in each group as values
snpGroupsDict = {} # key = snpGroup, value = array of snps in the group

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

# create isQTL_tissue feature vectors
isQTL_tissue_dict = {} # key = tissue feature name, value = isQTL dictionary for each tissue (key = snpGroup, value = isQTL)
for tissue in eQTL_dict:
	# for each snpGroup, determine whether any of its snps are or are associated with eQTLs
	isQTL_dict = {} # key = snpGroup, value = 1 or 0 (1 if any of the snps in the snpGroup is/is associated with an eQTL)
	num_isQTLs = 0

	for group in snpGroupsDict:
		snpsInEachGroup = snpGroupsDict[group]

		isQTL_values_inSnpGroup = []
		for i in range(len(snpsInEachGroup)):
			if (snpsInEachGroup[i] in eQTL_dict[tissue]) and (snpsInEachGroup[i] in gwas_isQTL_dict):
				isQTL_values_inSnpGroup.append(1)
			else:
				isQTL_values_inSnpGroup.append(0)

		if (1 in isQTL_values_inSnpGroup):
			isQTL_dict[group] = 1
			num_isQTLs += 1

		else:
			isQTL_dict[group] = 0

	print "Out of the total of", len(isQTL_dict), "GWAS snp groups in the model, there are", num_isQTLs, "snps that colocalize with eQTLs in", tissue

	tissue_feature = "isQTL_" + tissue
	isQTL_tissue_dict[tissue_feature] = isQTL_dict

# create output files
tab = "\t"
newline = "\n"

# create headerline
headerline = "snp" + tab
tissueList = isQTL_tissue_dict.keys()
tissueList.sort()
numTissueFeatures = len(tissueList)

for i in range(numTissueFeatures):
	if (i < (numTissueFeatures - 1)):
		headerline += tissueList[i] + tab
	else:
		headerline += tissueList[i] + newline

output = headerline

for snpGroup in snpGroupsDict:
	output += snpGroup + tab
	for i in range(numTissueFeatures):
		tissue = tissueList[i]
		if (i < (numTissueFeatures - 1)):
			output += str(isQTL_tissue_dict[tissue][snpGroup]) + tab
		else:
			output += str(isQTL_tissue_dict[tissue][snpGroup]) + newline

output_filename = "isQTL_byTissue.txt"
output_file = open(output_filename, 'w')
output_file.write(output)
output_file.close()