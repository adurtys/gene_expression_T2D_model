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

# store headers
headerLine = tissue_eqtl_file.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

for line in tissue_eqtl_file:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	eqtlInfo = []

	gene_id = columns[0]
	gene_chr = columns[2]
	var_tss_dist = columns[11]
	var_chr = columns[12]
	var_snp = columns[13]
	qval = columns[27]

	# only include if cis-eqtl and not on x-chromosome
	if (gene_chr != "X") and (int(gene_chr) == int(var_chr)):
		qvalList.append(float(qval))

		eqtlInfo.append(gene_id)
		eqtlInfo.append(int(gene_chr))
		eqtlInfo.append(var_tss_dist)
		eqtlInfo.append(int(var_chr))
		eqtlInfo.append(var_snp)
		eqtlInfo.append(float(qval))

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

print qtl_output_dict

# create output for output file
tab = "\t"
newline = "\n"

output = headers[0] + tab + headers[2] + tab + headers[11] + tab + headers[12] + tab + headers[27] + newline
for gene in qtl_output_dict:
	for i in range(6):
		if i < 5:
			output += str(qtl_output_dict[gene][i]) + tab
		else:
			output += str(qtl_output_dict[gene][i]) + newline

# create output file
outFilename = "sentinelQTL_test.txt"
outFile = open(outFilename, 'w')

outFile.write(output)

outFile.close()