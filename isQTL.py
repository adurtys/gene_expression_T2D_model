# Date Created: 16 November 2018
# Date Last Modified: 28 November 2018
# Execution: python isQTL.py [1] [2] [3]
# arvg[1] = matrix for eQTL snps that are in LD (r2 > 0.8) w/ GWAS snps (eQTL_inLD_wGWAS.txt)
# argv[2] = tissue information matrix for significant (FDR < 0.05) eQTL snps (signif_QTL_snps_allTisuses.txt)
# Description: creates matrix for whether each eQTL is in LD with a GWAS snp for every tissue

#!/usr/bin/env python
import sys

# read in command line arguments
eQTL_inLD_wGWAS_filename = sys.argv[1]
signif_eQTLs_byTissue_filename = sys.argv[2]

# read in eQTL snps in LD w/ GWAS snps
eQTL_inLD_wGWAS_file = open(eQTL_inLD_wGWAS_filename, 'r')

# skip headerline
eQTL_inLD_wGWAS_file.readline()

eQTLs_inLD_dictionary = {} # key = eQTL rsID, value = [eQTL, GWAS, GWAS snpGroup, chr]
for line in eQTL_inLD_wGWAS_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	eQTL = columns[1]
	gwas_snp = columns[2]
	gwas_snpGroup = columns[3]
	chromosome = columns[0]

	eQTLs_inLD_dictionary[eQTL] = [eQTL, gwas_snp, gwas_snpGroup, chromosome]

eQTL_inLD_wGWAS_file.close()

# read in signif_eQTLs by tissue
signif_eQTLs_byTissue_file = open(signif_eQTLs_byTissue_filename, 'r')

# skip headerline
signif_eQTLs_byTissue_file.readline()

eQTL_snplist = set()
eQTLs_byTissue_dictionary = {} # key = tissue, value = dictionary of eQTL snps and info for that tissue
for line in signif_eQTLs_byTissue_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	tissue = columns[0]
	eQTL = columns[6]
	eQTL_gene = columns[1]
	eQTL_chr = "chr" + columns[4]

	# add eqtl to snplist
	eQTL_snplist.add(eQTL)

	if tissue not in eQTLs_byTissue_dictionary:
		# create dictionary for that tissue
		eQTLs_byTissue_dictionary[tissue] = {} # key = eQTL, value = [eQTL, eQTL_chrNum, eQTL_gene]

	# add eQTL to tissue dictionary
	eQTLs_byTissue_dictionary[tissue][eQTL] = [eQTL, eQTL_chr, eQTL_gene]

signif_eQTLs_byTissue_file.close()

# create matrix for whether eQTL snp is in LD w/ a gwas snp for every tissue
isQTL_matrix = {} # key = eQTL, value = isQTL vector for each tissue (where isQTL = 1 if eQTL is in LD w/ a GWAS snp for that tissue)

tissueList = eQTLs_byTissue_dictionary.keys()
tissueList.sort()

# create dictionary of genes attached to all eQTLs that have isQTL = 1
isQTL_geneDict = {} # key = eQTL if isQTL == 1, value = gene attached to the eQTL that is in LD w/ a GWAS snp (assumes that only one gene is attached to each eQTL (TODO: CHECK IF THIS IS TRUE))

for eQTL in eQTL_snplist:
	isQTL_vector = []
	for tissue in tissueList:
		isQTL = 0
		if eQTL in eQTLs_byTissue_dictionary[tissue]:
			if eQTL in eQTLs_inLD_dictionary:
				isQTL = 1

				# obtain the gene attached to the eQTL
				eQTL_gene = eQTLs_byTissue_dictionary[tissue][eQTL][2]

				# store the gene attached to the eQTL
				if eQTL not in isQTL_geneDict:
					isQTL_geneDict[eQTL] = [eQTL_gene, tissue]
				else:
					# check that the stored gene is the same as the new gene
					if eQTL_gene != isQTL_geneDict[eQTL]:
						isQTL_geneDict[eQTL].append(eQTL_gene)
						isQTL_geneDict[eQTL].append(tissue)

		isQTL_vector.append(isQTL)

	isQTL_matrix[eQTL] = isQTL_vector

# create output files
tab = "\t"
newline = "\n"

isQTL_matrix_filename = "isQTL_matrix.txt"
isQTL_genes_filename = "isQTL_geneList.txt"

# create matrix
isQTL_matrix_file = open(isQTL_matrix_filename, 'w')

# create header line
header = "eQTL" + tab
for i in range(len(tissueList)):
	if (i < (len(tissueList) - 1)):
		header += tissueList[i] + tab
	else:
		header += tissueList[i] + newline

output = header

for eQTL in isQTL_matrix:
	output += eQTL + tab
	for i in range(len(isQTL_matrix[eQTL])):
		if (i < len(isQTL_matrix[eQTL]) - 1):
			output += str(isQTL_matrix[eQTL][i]) + tab
		else:
			output += str(isQTL_matrix[eQTL][i]) + newline

isQTL_matrix_file.write(output)
isQTL_matrix_file.close()

# create gene list
isQTL_genes_file = open(isQTL_genes_filename, 'w')

output = "eQTL" + tab + "ENSGID" + newline

for eQTL in isQTL_geneDict:
	output += eQTL + tab
	for i in range(len(isQTL_geneDict[eQTL])):
		if (i < len(isQTL_geneDict[eQTL]) - 1):
			output += isQTL_geneDict[eQTL][i] + tab
		else:
			output += isQTL_geneDict[eQTL][i] + newline

isQTL_genes_file.write(output)
isQTL_genes_file.close()