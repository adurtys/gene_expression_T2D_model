# Date Created: 16 November 2018
# Date Last Modified: 4 May 2019
# Execution: python isQTL.py [1] [2]
# argv[1] = matrix for eQTL snps that are in LD (r2 > 0.7) w/ GWAS snps (isLD.txt)
# argv[2] = grouped snps file (final_grouped_snps_all.txt)
# Description: creates ML table for whether each locus in the model is in LD with an eQTL snp

#!/usr/bin/env python
import sys

# read in command line arguments
eQTL_inLD_wGWAS_filename = sys.argv[1]

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
groupedSnpsFilename = sys.argv[2]

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

# for each snpGroup, determine whether any of its snps are or are associated with eQTLs
isQTL_dict = {} # key = snpGroup, value = 1 or 0 (1 if any of the snps in the snpGroup is/is associated with an eQTL)

num_isQTLs = 0
for group in snpGroupsDict:
	snpsInEachGroup = snpGroupsDict[group]
	
	isQTL_values_inSnpGroup = []
	for i in range(len(snpsInEachGroup)):
		if snpsInEachGroup[i] in gwas_isQTL_dict:
			isQTL_values_inSnpGroup.append(1)
		else:
			isQTL_values_inSnpGroup.append(0)

	if (1 in isQTL_values_inSnpGroup):
		isQTL_dict[group] = 1
		num_isQTLs += 1
	else:
		isQTL_dict[group] = 0

print "Out of the total of", len(isQTL_dict), "GWAS snp groups in the model, there are", num_isQTLs, "snps that colocalize with eQTLs."

# create output files
tab = "\t"
newline = "\n"

output = "snp" + tab + "isQTL" + newline

for snpGroup in isQTL_dict:
	output += snpGroup + tab + str(isQTL_dict[snpGroup]) + newline

output_filename = "isQTL.txt"
output_file = open(output_filename, 'w')
output_file.write(output)
output_file.close()



## then create separate features for each of the 53 tissues

## then create the ML table



# eQTL_snplist = set()
# for line in signif_eQTLs_byTissue_file:


# 	tissue = columns[0]
# 	eQTL = columns[6]
# 	eQTL_gene = columns[1]
# 	eQTL_chr = "chr" + columns[4]

# 	# add eqtl to snplist
# 	eQTL_snplist.add(eQTL)

# 	if tissue not in eQTLs_byTissue_dictionary:
# 		# create dictionary for that tissue
# 		eQTLs_byTissue_dictionary[tissue] = {} # key = eQTL, value = [eQTL, eQTL_chrNum, eQTL_gene]

# 	# add eQTL to tissue dictionary
# 	eQTLs_byTissue_dictionary[tissue][eQTL] = [eQTL, eQTL_chr, eQTL_gene]

# signif_eQTLs_byTissue_file.close()

# # create matrix for whether eQTL snp is in LD w/ a gwas snp for every tissue
# isQTL_matrix = {} # key = eQTL, value = isQTL vector for each tissue (where isQTL = 1 if eQTL is in LD w/ a GWAS snp for that tissue)

# tissueList = eQTLs_byTissue_dictionary.keys()
# tissueList.sort()

# # create dictionary of genes attached to all eQTLs that have isQTL = 1
# isQTL_geneDict = {} # key = eQTL if isQTL == 1, value = gene attached to the eQTL that is in LD w/ a GWAS snp (assumes that only one gene is attached to each eQTL (TODO: CHECK IF THIS IS TRUE))

# for eQTL in eQTL_snplist:
# 	isQTL_vector = []
# 	for tissue in tissueList:
# 		isQTL = 0
# 		if eQTL in eQTLs_byTissue_dictionary[tissue]:
# 			if eQTL in eQTLs_inLD_dictionary:
# 				isQTL = 1

# 				# obtain the gene attached to the eQTL
# 				eQTL_gene = eQTLs_byTissue_dictionary[tissue][eQTL][2]

# 				# store the gene attached to the eQTL
# 				if eQTL not in isQTL_geneDict:
# 					isQTL_geneDict[eQTL] = [eQTL_gene, tissue]
# 				else:
# 					# check that the stored gene is the same as the new gene
# 					if eQTL_gene != isQTL_geneDict[eQTL][0]:
# 						isQTL_geneDict[eQTL].append(eQTL_gene)
# 						isQTL_geneDict[eQTL].append(tissue)

# 		isQTL_vector.append(isQTL)

# 	isQTL_matrix[eQTL] = isQTL_vector

# # create output files
# tab = "\t"
# newline = "\n"

# isQTL_matrix_filename = "isQTL_matrix.txt"
# isQTL_genes_filename = "isQTL_geneList.txt"

# # create matrix
# isQTL_matrix_file = open(isQTL_matrix_filename, 'w')

# # create header line
# header = "eQTL" + tab
# for i in range(len(tissueList)):
# 	if (i < (len(tissueList) - 1)):
# 		header += tissueList[i] + tab
# 	else:
# 		header += tissueList[i] + newline

# output = header

# for eQTL in isQTL_matrix:
# 	output += eQTL + tab
# 	for i in range(len(isQTL_matrix[eQTL])):
# 		if (i < len(isQTL_matrix[eQTL]) - 1):
# 			output += str(isQTL_matrix[eQTL][i]) + tab
# 		else:
# 			output += str(isQTL_matrix[eQTL][i]) + newline

# isQTL_matrix_file.write(output)
# isQTL_matrix_file.close()

# # create gene list
# isQTL_genes_file = open(isQTL_genes_filename, 'w')

# output = "eQTL" + tab + "ENSGID" + tab + "TISSUE" + newline

# for eQTL in isQTL_geneDict:
# 	output += eQTL + tab
# 	for i in range(len(isQTL_geneDict[eQTL])):
# 		if (i < len(isQTL_geneDict[eQTL]) - 1):
# 			output += isQTL_geneDict[eQTL][i] + tab
# 		else:
# 			output += isQTL_geneDict[eQTL][i] + newline

# isQTL_genes_file.write(output)
# isQTL_genes_file.close()