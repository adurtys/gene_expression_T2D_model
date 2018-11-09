# Date Created: 2 November 2018
# Date Last Modified: 9 November 2018
# Execution: python eQTL_gwas_ld_snps.py [1] [2]
# argv[1] = path to directory containing LD files for all chromosomes (/project/voight_ML/adurtys/eqtl_feature/isLD)
# argv[2] = list of eQTL snps
# Description: determines which eQTL snps are in LD with GWAS snps of interest
# Run Time: 4 sec

#!/usr/bin/env python
import sys, os

# read in command-line arguments
ld_directory = sys.argv[1]
eQTL_snp_filename = sys.argv[2]

ld_filenames = os.listdir(ld_directory)

all_ld_snps_dict = {} # key = chromosome, value = dictionary of eQTL/GWAS snps in LD with each other for that chromosome

for filename in ld_filenames:
	# store choromosome number
	chromosome = filename.split('.')[0].strip("isLD_")
	filepath = ld_directory + "/" + filename

	file = open(filepath, 'r')

	ld_snps_dict = {} # key = snpA, value = [snpA, snpB, r2]
	ld_snp_pairs = {}

	for line in file:
		line = line.rstrip('\r\n')
		columns = line.split()

		snpA = columns[2]
		snpB = columns[5]
		r2 = float(columns[6])

		if (snpA != snpB):
			if snpA not in ld_snp_pairs.keys():
				# check if it is the values (stored as snpB)
				if snpA in ld_snp_pairs.values():
					# snpA is stored as snpB in the pairs dictionary --> get snpA
					origSnpA = ld_snp_pairs.keys()[ld_snp_pairs.values().index(snpA)]
					# check that r2 is the same
					if (ld_snps_dict[origSnpA][2] != r2):
						print "ERROR: r2 values are different for", snpA

				else: #store snpA
					ld_snps_dict[snpA] = snpB
					ld_snps_dict[snpA] = [snpA, snpB, r2, chromosome]

	all_ld_snps_dict[chromosome] = ld_snps_dict

	file.close()
	
eQTL_snp_file = open(eQTL_snp_filename, 'r')
eQTL_snpList = []
for line in eQTL_snp_file:
	eQTL_snp = line.rstrip('\r\n')
	eQTL_snpList.append(eQTL_snp)
eQTL_snp_file.close()

tab = "\t"
newline = "\n"

output = ""

for snp in eQTL_snpList:
	for chromosome in all_ld_snps_dict:
		if (snp in all_ld_snps_dict[chromosome]):
			for i in range(len(all_ld_snps_dict[chromosome][snp])):
				if i < (len(all_ld_snps_dict[chromosome][snp]) - 1):
					output += str(all_ld_snps_dict[chromosome][snp][i]) + tab
				else:
					output += str(all_ld_snps_dict[chromosome][snp][i]) + newline

outFilename = eQTL_inLD_wGWAS.txt
outFile = open(outFilename, 'w')
outFile.write(output)
outFile.close()