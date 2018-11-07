# Date Created: 2 November 2018
# Date Last Modified: 2 November 2018
# Execution: python eQTL_gwas_ld_snps.py [1] [2]
# argv[1] = path to directory containing LD files for all chromosomes (/project/voight_ML/adurtys/eqtl_feature/plink_output/ld_files)
# argv[2] = eQTL_snplist_filename (uniq_eQTL_snp_rsIDs.txt)
# Description: determines which eQTL snps are in LD with GWAS snps of interest
# Run Time: 

#!/usr/bin/env python
import sys, os

# read in command-line arguments
ld_directory = sys.argv[1]
eQTL_snplist_filename = sys.argv[2]

ld_filenames = os.listdir(ld_directory)

all_ld_snps_dict = {} # key = chromosome, value = dictionary of eQTL/GWAS snps in LD with each other for that chromosome

for filename in ld_filenames:
	# store choromosome number
	chromosome = filename.split('.')[0].strip("LD_")
	filepath = ld_directory + "/" + filename

	file = open(filepath, 'r')

	file.readline() # remove header line

	ld_snps_dict = {} # key = snpA, value = [snpA, snpB, r2]

	for line in file:
		line = line.rstrip('\r\n')
		columns = line.split('\t')

		snpA = columns[2]
		snpB = columns[5]
		r2 = float(columns[6])

		if (r2 == 1) and (snpA != snpB):
			print "ERROR: r2 = 1, but snps A and B are not the same! snpA:", snpA, "snpB:", snpB
		elif (r2 != 1):
			if (snpA not in ld_snps_dict) and (snpA not in ld_snps_dict.values()[1]):
				ld_snps_dict[snpA] = 