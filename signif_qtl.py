# Date Created: 26 October 2018
# Date Last Modified: 18 February 2019
# Execution: python signif_qtl.py [1] [2]
# argv[1] = path to directory containing GTEx tissue eQTL files (/project/voight_datasets/GTEx_V6p/egenes)
# argv[2] = threshold for q-value significance for variant-gene pair (0.05)
# Description: creates a matrix of qtls that are significant (q-value less than threshold) for all tissues
# Run time: ~20 sec

#!/usr/bin/env python
import sys, os

# read in command-line arguments
tissueQTL_Directory = sys.argv[1]
threshold = float(sys.argv[2])

tissueQTL_filenames = os.listdir(tissueQTL_Directory)

# key = tissue type, value = dictionary of info QTLs in LD with a GWAS snp on the same chromosome and within a certain SPECIFIED distance (TODO)
tissue_signifQTL_Dict = {} 

for filename in tissueQTL_filenames:
	# store tissue name
	tissue = filename.rstrip("Analysis.v6p.egenes.txt").rstrip('_')
	filepath = tissueQTL_Directory + "/" + filename
	
	file = open(filepath, 'r')

	# store the header line once (for first file in the directory)
	if (filename == tissueQTL_filenames[0]):
		headerLine = file.readline()
	else:
		file.readline()

	eqtlDict = {}

	for line in file:
		line = line.rstrip('\r\n')
		columns = line.split('\t')

		eqtlInfo = []

		gene_id = columns[0]
		gene_chr = columns[2]
		var_tss_dist = columns[11]
		var_chr = columns[12]
		var_snp = columns[13]
		rsID = columns[16]
		qval = float(columns[27])

		# only include if cis-eqtl and not on x-chromosome
		if (gene_chr != "X") and (int(gene_chr) == int(var_chr)):
			if (qval < threshold):
				eqtlInfo.extend([tissue, gene_id, int(gene_chr), int(var_tss_dist), int(var_chr), var_snp, rsID, qval])

		eqtlDict[gene_id] = eqtlInfo

	tissue_signifQTL_Dict[tissue] = eqtlDict

	file.close()

# create output directory
outputDirectoryPath = "./signif_QTL_snps_directory"
os.mkdir(outputDirectoryPath)

# change directory to output directory
os.chdir(outputDirectoryPath)

tab = "\t"
newline = "\n"

# create new header line
headers = headerLine.split('\t')

newHeaderLine = "Tissue" + tab + headers[0] + tab + headers[2] + tab + headers[11] + tab + headers[12] + tab + headers[13] + tab + headers[16] + tab + headers[27] + newline

# create output file for all tissues (combined)
combinedOutputFilename = "signif_QTL_snps_allTisuses.txt"
combinedOutputFile = open(combinedOutputFilename, 'w')

combinedOutput = newHeaderLine
numTotalQTL = 0
for tissue in tissue_signifQTL_Dict:
	# create output file for each tissue
	outputFilename = "signif_QTL_snps_" + tissue
	outputFile = open(outputFilename, 'w')
	
	output = newHeaderLine
	numTissueQTL = 0
	for gene_id in tissue_signifQTL_Dict[tissue]:
		for i in range(len(tissue_signifQTL_Dict[tissue][gene_id])):
			if (i < len(tissue_signifQTL_Dict[tissue][gene_id]) - 1):
				output += str(tissue_signifQTL_Dict[tissue][gene_id][i]) + tab
				combinedOutput += str(tissue_signifQTL_Dict[tissue][gene_id][i]) + tab
			else: # add newline at end of each gene
				output += str(tissue_signifQTL_Dict[tissue][gene_id][i]) + newline
				numTissueQTL += 1
				combinedOutput += str(tissue_signifQTL_Dict[tissue][gene_id][i]) + newline
				numTotalQTL += 1

	outputFile.write(output)
	outputFile.close()

	print tissue, "has", numTissueQTL, "cis-eqtls associated for FDR <", threshold

combinedOutputFile.write(combinedOutput)
combinedOutputFile.close()

print "There are a total of", numTotalQTL, "cis-eqtls across all", len(tissue_signifQTL_Dict), "tissues."