# Date Created: 2 November 2018
# Date Last Modified: 29 April 2019
# Execution: python isLD.py [1] [2] [3]
# argv[1] = linked snps (A & B) for all chromosomes (./plink/gwas_eQTL_LD)
# argv[2] = list of eQTL snps (./uniq_eQTL_snp_rsIDs.txt)
# argv[3] = list of GWAS snps (./uniq_GWAS_snp_rsIDs.txt)
# Description: determines which eQTL snps are in LD with GWAS snps of interest
# Run Time: 2 sec

#!/usr/bin/env python
import sys, os

# read in command-line arguments
ld_filename = sys.argv[1]
eQTL_snp_filename = sys.argv[2]
gwas_snp_filename = sys.argv[3]

# read in list of eQTL snps
eQTL_snp_file = open(eQTL_snp_filename, 'r')
eQTL_rsID_list = set()
for line in eQTL_snp_file:
	eQTL_rsID = line.rstrip('\r\n')
	eQTL_rsID_list.add(eQTL_rsID)
eQTL_snp_file.close()
print "Finished reading in eQTL_snp_file. This file has", len(eQTL_rsID_list), "snps."

# read in list of GWAS snps
gwas_snp_file = open(gwas_snp_filename, 'r')
gwas_rsID_list = set()
for line in gwas_snp_file:
	gwas_rsID = line.rstrip('\r\n')
	gwas_rsID_list.add(gwas_rsID)
gwas_snp_file.close()
print "Finished reading in gwas_snp_file. This file has", len(gwas_rsID_list), "snps."

# create output file for every eQTL snp that is in LD with a GWAS snp
tab = "\t"
newline = "\n"

output = "CHR_A_BP_A_" + tab + "SNP_A" + tab + "CHR_B_BP_B_" + tab + "SNP_B" + newline

# read in ld file
ld_file = open(ld_filename,'r')

# dictionary to store eQTLs that are in LD with GWAS snps --> key = snpA (eQTL snp), value = [eQTL snp, snpB_position, snpB_rsID]
isLD_dict = {}
isLD_multipleLinkedSnps_dict = {} # key = rsID for eQTL snp linked to multiple gwas snps, value = [GWAS snps (position and rsID) linked to the eQTL snp]

ld_file.readline()

numColocalizedSnps = 0
for line in ld_file:
	line = line.rstrip('\r\n')
	columns = line.split()

	snpA_position = columns[0]
	snpA_rsID = columns[1]
	snpB_position = columns[2]
	snpB_rsID = columns[3]

	if (snpA_rsID in gwas_rsID_list) and (snpB_rsID in eQTL_rsID_list):
		numColocalizedSnps += 1

		output += snpA_position + tab + snpA_rsID + tab + snpB_position + tab + snpB_rsID + newline

outFilename = "isLD_2.txt"
outFile = open(outFilename, 'w')
outFile.write(output)
outFile.close()