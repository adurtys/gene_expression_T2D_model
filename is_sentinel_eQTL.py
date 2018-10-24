# Date Created: 24 October 2018
# Date Last Modified: 24 October 2018
# Execution: python is_sentinel_eQTL.py 
# argv1: tissue_eqtl.txt
# Description: for the given tissue file, outputs the 5 qtls with the smallest false discovery rate (q-val) and the snps attached to them
# Run time: 

#!/usr/bin/env python
import sys

# read in command line argument
tissue_eqtl_filename = sys.argv[1]

# read in tissue_eqtl file
tissue_eqtl_file = open(tissue_eqtl_filename, 'r')

# create data structures to store info from the file
eqtlDict = {}
qvalList = []

for line in tissue_eqtl_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	eqtlInfo = []

	gene_id = columns[0]
	gene_chr = int(columns[1])
	var_tss_dist = columns[11]
	var_chr = int(columns[12])
	qval = columns[27]

	# only include if cis-eqtl
	if (gene_chr == var_chr):
		qvalList.append(float(qval))

	eqtlInfo.append(gene_id)
	eqtlInfo.append(gene_chr)
	eqtlInfo.append(var_tss_dist)
	eqtlInfo.append(var_chr)
	eqtlInfo.append(qval)

	eqtlDict[gene_id] = eqtlInfo

tissue_eqtl_file.close()

qvalList.sort()

# determine five smallest q-values
smallest_qval = []
for i in range(5):
	smallest_qval.append(qvalList[i])

# determine the genes attached to the smallest q-values
qtl_output_dict = {}
for gene in eqtlDict:
	if eqtlDict[gene][4] in smallest_qval:
		qtl_output_dict[gene] = eqtlDict[gene]

# create output for output file
tab = "\t"
newline = "\t"

output = ""
for gene in qtl_output_dict:
	for i in range(5):
		if i < 4:
			output += qtl_output_dict[gene][i] + tab
		else:
			output += qtl_output_dict[gene][i] + newline

# create output file
outFilename = "sentinelQTL_test.txt"
outFile = open(outFilename, 'w')

outFile.write(output)

outFile.close()