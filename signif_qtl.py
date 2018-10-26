# Date Created: 26 October 2018
# Date Last Modified: 26 October 2018
# Execution: python signif_qtl.py [1] [2]
# argv[1] = path to directory containing GTEx tissue eQTL files (/project/voight_datasets/GTEx_V6p/egenes)
# argv[2] = threshold for q-value significance for variant-gene pair (0.05)
# Description: creates a matrix of qtls that are significant (q-value less than threshold) for all tissues and a "collapsed tissue"
# Run time: 

import sys, os

# read in command-line arguments
tissueQTL_Directory = sys.argv[1]
threshold = float(sys.argv[2])

tissueQTL_filenames = os.listdir(tissueQTL_Directory)

# key = tissue type, value = dictionary of info QTLs in LD with a GWAS snp on the same chromosome and within a certain SPECIFIED distance (TODO)
tissue_signifQTL_Dict = {} 

for filename in tissueQTL_files:
	# store tissue name
	tissue = filename.rstrip("Analysis.v6p.egenes.txt").rstrip('_')
	
	file = open(filename, 'r')

	# store the header line once (for first file in the directory)
	if (filename = tissueQTL_files[0]):
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
				eqtlInfo.extend([gene_id, int(gene_chr), int(var_tss_dist), int(var_chr), var_snp, rsID, qval])

		eqtlDict[gene_id] = eqtlInfo

	tissue_signifQTL_Dict[tissue] = eqtlDict

	file.close()

# create output file
outputFilename = "signif_QTL_snps_allTisuses.txt"
outputFile = open(outputFilename, 'w')

tab = "\t"
newline = "\n"

# create new header line
newHeaderLine = "Tissue" + tab + headerLine

output = newHeaderLine
for tissue in tissue_signifQTL_Dict:
	for gene_id in tissue_signifQTL_Dict[tissue]:
		output += tissue + tab
		for i in range(len(tissue_signifQTL_Dict[tissue][gene_id])):
			if (i < len(tissue_signifQTL_Dict[tissue][gene_id]) - 1):
				output += str(tissue_signifQTL_Dict[tissue][gene_id][i]) + tab
			else: # add newline at end of each gene
				output += str(tissue_signifQTL_Dict[tissue][i]) + newline

outputFile.write(output)
outputFile.close()