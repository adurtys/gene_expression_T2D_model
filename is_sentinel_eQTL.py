# Date Created: 24 October 2018
# Date Last Modified: 24 October 2018
# Execution: python is_sentinel_eQTL.py 
# argv1: tissue_eqtl.txt
# Description: for the given tissue file, outputs the 5 qtls with the smallest p-val and the snps attached to them
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
pvalList = []

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
	rsID = columns[16]
	pval = columns[22]
	qval = columns[27]

	# only include if cis-eqtl and not on x-chromosome
	if (gene_chr != "X") and (int(gene_chr) == int(var_chr)):
		qvalList.append(float(qval))
		pvalList.append(float(pval))

		eqtlInfo.append(gene_id)
		eqtlInfo.append(int(gene_chr))
		eqtlInfo.append(var_tss_dist)
		eqtlInfo.append(int(var_chr))
		eqtlInfo.append(var_snp)
		eqtlInfo.append(rsID)
		eqtlInfo.append(float(pval))
		eqtlInfo.append(float(qval))

		eqtlDict[gene_id] = eqtlInfo

tissue_eqtl_file.close()

qvalList.sort()
pvalList.sort()

# determine five smallest q-values
smallest_qval = []
for i in range(5):
	smallest_qval.append(qvalList[i])

# determine the genes attached to the smallest q-values
qval_qtl_output_dict = {}
for gene in eqtlDict:
	if eqtlDict[gene][7] in smallest_qval:
		qval_qtl_output_dict[gene] = eqtlDict[gene]

print qval_qtl_output_dict

# create output for output file
tab = "\t"
newline = "\n"

qval_output = headers[0] + tab + headers[2] + tab + headers[11] + tab + headers[12] + tab + headers[13] + tab + headers[16] + tab + headers[22] + tab + headers[27] + newline
for gene in qval_qtl_output_dict:
	for i in range(8):
		if i < 7:
			qval_output += str(qval_qtl_output_dict[gene][i]) + tab
		else:
			qval_output += str(qval_qtl_output_dict[gene][i]) + newline

# create output file
qval_outFilename = "sentinelQTL_qval_test.txt"
qval_outFile = open(qval_outFilename, 'w')

qval_outFile.write(qval_output)

qval_outFile.close()

# determine five smallest p-values
smallest_pval = []
for i in range(5):
	smallest_pval.append(pvalList[i])

# determine the genes attached to the smallest q-values
pval_qtl_output_dict = {}
for gene in eqtlDict:
	if eqtlDict[gene][6] in smallest_pval:
		pval_qtl_output_dict[gene] = eqtlDict[gene]

print pval_qtl_output_dict

pval_output = headers[0] + tab + headers[2] + tab + headers[11] + tab + headers[12] + tab + headers[13] + tab + headers[22] + tab + headers[27] + newline
for gene in pval_qtl_output_dict:
	for i in range(8):
		if i < 7:
			pval_output += str(pval_qtl_output_dict[gene][i]) + tab
		else:
			pval_output += str(pval_qtl_output_dict[gene][i]) + newline

# create output file
pval_outFilename = "sentinelQTL_pval_test.txt"
pval_outFile = open(pval_outFilename, 'w')

pval_outFile.write(pval_output)

pval_outFile.close()