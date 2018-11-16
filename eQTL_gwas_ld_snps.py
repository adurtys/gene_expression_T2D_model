# Date Created: 2 November 2018
# Date Last Modified: 16 November 2018
# Execution: python eQTL_gwas_ld_snps.py [1] [2] [3] [4]
# argv[1] = path to directory containing LD files for all chromosomes (/project/voight_ML/adurtys/eqtl_feature/plink_output/ld_files)
# argv[2] = list of eQTL snps (uniq_eQTL_snp_rsIDs.txt)
# argv[3] = list of GWAS snps (all_GREGOR_snplist_rsIDs.txt)
# argv[4] = threshold for LD significance (r2 > 0.8)
# Description: determines which eQTL snps are in LD with GWAS snps of interest
# Run Time: 4 sec

#!/usr/bin/env python
import sys, os

# read in command-line arguments
ld_directory = sys.argv[1]
eQTL_snp_filename = sys.argv[2]
gwas_snp_filename = sys.argv[3]
threshold = float(sys.argv[4])

# read in list of eQTL snps
eQTL_snp_file = open(eQTL_snp_filename, 'r')
eQTL_rsID_list = set()
for line in eQTL_snp_file:
	eQTL_rsID = line.rstrip('\r\n')
	
	eQTL_rsID_list.add(eQTL_rsID)
eQTL_snp_file.close()

# read in list of GWAS snps
gwas_snp_file = open(gwas_snp_filename, 'r')
gwas_snp_dict = {} # key = rsID, value = [chrNum, chr#:group]
for line in gwas_snp_file:
	line = line.rstrip('\r\n')
	columns = line.split()

	gwas_rsID = columns[0]
	snpGroup = columns[1]
	chrNum = int(snpGroup.split(':')[0].strip("chr"))

	gwas_snp_dict[gwas_rsID] = [chrNum, snpGroup]
gwas_snp_file.close()

ld_filenames = os.listdir(ld_directory)

all_eQTL_inLD = {} # key = chromosome, value = dictionary of eQTL/GWAS snps in LD with each other for that chromosome

for filename in ld_filenames:
	# store choromosome number
	chromosome = filename.split('.')[0].strip("LD_")
	filepath = ld_directory + "/" + filename

	file = open(filepath, 'r')
	file.readline()

	eQTL_inLD_wGWAS = {} # key = rsID of eQTL snp, value = GWAS_rsID, r2

	for line in file:
		line = line.rstrip('\r\n')
		columns = line.split()

		snpA = columns[2]
		snpB = columns[5]
		r2 = float(columns[6])

		if (snpA in eQTL_rsID_list) and (snpB in gwas_snp_dict) and (r2 > threshold):
			# snpA is an eQTl that is in LD with a gwas snp
			eQTL_inLD_wGWAS[snpA] = [snpB, gwas_snp_dict[snpB][0], gwas_snp_dict[snpB][1], r2]

	all_eQTL_inLD[chromosome] = eQTL_inLD_wGWAS

	file.close()

# create output file for every eQTL snp (by chromosome) that is in LD with a GWAS snp
tab = "\t"
newline = "\n"

output = ""

for chromosome in all_eQTL_inLD:
	for eQTL_rsID in all_eQTL_inLD[chromosome]:
		output += chromosome + tab + eQTL_rsID + tab
		for i in range(len(all_eQTL_inLD[chromosome][eQTL_rsID])):
			if i < (len(all_eQTL_inLD[chromosome][eQTL_rsID]) - 1):
				output += all_eQTL_inLD[chromosome][eQTL_rsID][i] + tab
			else:
				output += all_eQTL_inLD[chromosome][eQTL_rsID] + newline

outFilename = "eQTL_inLD_wGWAS.txt"
outFile = open(outFilename, 'w')
outFile.write(output)
outFile.close()
