# Date Created: 3 December 2018
# Date Last Modified: 3 December 2018
# Execution: python isNG.py [1] [2]
# [1]: final_grouped_snps_all.txt (from /project/voight_ML/lorenzk/V2/Public/models/snplist)
# [2]: all_GREGOR_snplist_rsIDs.txt (from /project/voight_ML/lorenzk/V2/Public/models/snplist)
# [3]: nearestGenes_centroidSnps_final_grouped_snps_all.txt (from rm_willer_diagram)
# [4]: isQTL_geneList.txt
# [5]: eQTL_tissues.txt
# Description: creates matrix for whether the gene attached to each eQTL snp that was in LD with a GWAS snp for every tissue is one of the nearest genes to that snp

#!/usr/bin/env python
import sys

# read in command line arguments
grouped_snp_filename = sys.argv[1]
snp_rsID_filename = sys.argv[2]
geneExp_ng_filename = sys.argv[3]
eQTL_geneList_filename = sys.argv[4]
eQTL_tissue_filename = sys.argv[5]

# create dictionary to store data for snp
snp_info_dict = {} # key = snp (chr:#), value = [snpGroup (Chr#_Group_#), rsID]

# store snp and its group
grouped_snp_file = open(grouped_snp_filename, 'r')
for line in grouped_snp_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snpGroup = columns[0]
	snp = columns[1]

	snp_info_dict[snp] = snpGroup
grouped_snp_file.close()

# add rsID to snp info dictionary 
snp_rsID_file = open(snp_rsID_filename, 'r')
for line in snp_rsID_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	rsID = columns[0]
	snp = columns[1]

	if snp in snp_info_dict:
		snp_info_dict[snp].append(rsID)

snp_rsID_file.close()

# store nearest genes from gene expression lookup results for each rsID
ng_dict = {} # key = snp, value = [firstNG, secondNG, thirdNG]

geneExp_ng_file = open(geneExp_ng_filename, 'r')
for line in geneExp_ng_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snpGroup = columns[0]
	firstNG = columns[1]
	secondNG = columns[2]
	thirdNG = columns[3]

	ng_dict[snpGroup] = [firstNG, secondNG, thirdNG]
geneExp_ng_file.close()

# create rsID dict to convert from rsID to snpGroup
rsID_to_snpGroup_dict = {}
for snp in snp_info_dict:
	rsID = snp_info_dict[snp][1]
	snpGroup = snp_info_dict[snp][0]

	rsID_to_snpGroup_dict[rsID] = snpGroup

# read in isQTL gene list
isQTL_tissue_dict = {} # key = groupedSnp, value = {tissue:gene}

eQTL_geneList_file = open(eQTL_geneList_filename, 'r')
for line in isQTL_geneList:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	rsID = columns[0]

	# obtain snpGroup corresponding to the rsID
	snpGroup = ""
	if rsID in rsID_to_snpGroup_dict:
		snpGroup = rsID_to_snpGroup_dict[rsID]
	else:
 		print "ERROR (isNG.py line 89):", rsID, "in eQTL_expMatrix_file, but does not contain a corresponding snpGroup"

	numCol = len(columns)

	# store tissue(s) and corresponding gene(s) in tissue dictionary
	for i in range(1, numCol): # ignore the fist column, which contains the rsID
		if (i % 2) != 0:
			gene = columns[i].split('.')[0]
			tissue = columns[i + 1]

			# update isQTL tissue dictionary
			if snpGroup not in isQTL_tissue_dict:
				isQTL_tissue_dict[snpGroup] = {tissue : gene}
			else:
				isQTL_tissue_dict[snpGroup][tissue] = gene
eQTL_geneList_file.close()

tissueList = []

# read in tissues
eQTL_tissue_file = open(eQTL_tissue_filename, 'r')
for line in eQTL_tissue_file:
	tissue = line.rstrip('\r\n')
	tissueList.append(tissue)
eQTL_tissue_file.close()

# create isNG matrix
isNG_matrix = {} # key = snpGroup, value = vector for whether, for each tissue, the snp is in LD with a GWAS snp, and the corresponding gene is one of the nearest genes

for snpGroup in isQTL_tissue_dict:
	isNG_vector = []

	for tissue in tissueList:
		isNG = 0
		if tissue in isQTL_tissue_dict[snpGroup]:
			geneToCompare = isQTL_tissue_dict[snp][tissue]
			nearestGenes = ng_dict[snpGroup]

			if geneToCompare in nearestGenes:
				isNG = 1
		isNG_vector.append(isNG)

	isNG_matrix[snpGroup] = isNG_vector

# create isNG output matrix
tab = "\t"
newline = "\n"

isNG_matrix_filename = "isNG_matrix.txt"
isNG_matrix_file = open(isNG_matrix_filename, 'w')

# create header line
header = ""

for i in range(len(tissueList)):
	if i < (len(tissueList) - 1):
		header += tissueList[i] + tab
	else:
		header += tissueList[i] + newline

output = header

for snpGroup in isNG_matrix:
	output += snpGroup + tab
	for i in range(len(isNG_matrix[snpGroup])):
		if i < (len(isNG_matrix[snpGroup]) - 1):
			output += str(isNG_matrix[snpGroup][i]) + tab
		else:
			output += str(isNG_matrix[snpGroup][i]) + newline

isNG_matrix_file.write(output)
isNG_matrix_file.close()