# Date Created: 14 November 2018
# Date Last Modified: 14 November 2018
# Execution: python isLD_ML_table.py [1] [2] [3] [4] [5]
# argv[1] = eQTL_inLD_wGWAS.txt (78916)
# argv[2] = GWAS_snps_inLD.txt (8890)
# argv[3] = all_GREGOR_snplist_rsIDs.txt (17797)
# argv[4] = signif_QTL_snps_allTissues.txt (194313)
# argv[5] = threshold for significant LD between eQTL snp and GWAS snp (= 0.8)
# Description: creates ML table --> for every gwas snp group, determines whether it is in LD with a QTL for every tissue 
# Runtime: 

#!/usr/bin/env python
import sys

# read in command-line arguments
eQTL_inLD_wGWAS_filename = sys.argv[1]
GWAS_snps_inLD_filename = sys.argv[2]
all_GWAS_snps_filename = sys.argv[3]
eQTLs_byTissue_filename = sys.argv[4]
threshold = float(sys.argv[5])

gwasInLD_snpList = set() # len = total number of GWAS in LD (not necessarily significant, but associated from plink output)
# read in GWAS snps that are in LD
GWAS_snps_inLD_file = open(GWAS_snps_inLD_filename, 'r')
for line in GWAS_snps_inLD_file:
	snp = line.rstrip('\r\n')
	gwasInLD_snpList.add(snp)
GWAS_snps_inLD_file.close()

# store info for eQTLs in LD w/ GWAS (if r2 > threshold)
eQTLs_inLD_dict = {}

eQTL_inLD_wGWAS_file = open(eQTL_inLD_wGWAS_filename, 'r')
for line in eQTL_inLD_wGWAS_file:
	line = line.rstrip('\r\n')
	columns = line.split()

	eQTL_snp = columns[0]
	GWAS_snp = columns[1]
	r2 = float(columns[2])
	chromosome = columns[3]

	if GWAS_snp in gwas_snpList:
		if r2 > threshold:
			eQTLs_inLD_dict[eQTL_snp] = [eQTL_snp, GWAS_snp, r2, chromosome]
eQTL_inLD_wGWAS_file.close()

print "Out of the" len(gwasInLD_snpList), "GWAS snps associated with eQTL snps,", len(eQTLs_inLD_dict), "are significantly in LD, with an r2 >", threshold

# read in all GWAS snp groups being used in model
all_GWAS_snps_file = open(all_GWAS_snps_filename, 'r')
GWAS_snp_dict = {} # key = snpGroup, value = rsID for snpGroup
for line in all_GWAS_snps_file:
	line = line.rstrip('\r\n')
	columns = line.split()

	rsID = columns[0]
	snpGroup = columns[1]

	GWAS_snp_dict[snpGroup] = rsID
all_GWAS_snps_file.close()

# read in eQTL tissue information
eQTLs_byTissue_file = open(eQTLs_byTissue_filename, 'r')
eQTLs_byTissue_matrix = {} # key = tissue, value = dictionary of eQTL info for every rsID for that tissue
for line in eQTLs_byTissue_file::
	line = line.rstrip('\r\n')
	columns = line.split()

	tissue = columns[0]
	rsID = columns[6]

	if tissue not in eQTLs_byTissue_matrix.keys():
		# create dictionary for the tissue
		eQTLs_byTissue_matrix[tissue] = {} # key = rsID, value = list of gene eqtl info
		eQTLs_byTissue_matrix[tissue][rsID] = []
		
		for i in len(columns):
			eQTLs_byTissue_matrix[tissue][rsID].append(columns[i])
	
	else: # tissue already has a corresponding dictionary
		eQTLs_byTissue_matrix[tissue][rsID] = []

		for i in len(columns):
			eQTLs_byTissue_matrix[tissue][rsID].append(columns[i])

eQTLs_byTissue_file.close()

# create vector dictionary for whether GWAS 


