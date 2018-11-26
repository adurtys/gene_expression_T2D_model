#=======================================================================================================
# README for adding gene expression features to the ML table
#=======================================================================================================

This file describes the inputs required for files that should be executed in order to append gene expression features to the T2D model. 

#-------------------------------------------------------------------------------------------------------
# to obtain centroid snps
#-------------------------------------------------------------------------------------------------------

Execution: python centroidPositions_willerDiagram.py [1]
	[1] = filename for file containing grouped snps

Output: 
	centroidSnps_final_grouped_snps_all.txt

#-------------------------------------------------------------------------------------------------------
# to conduct gene expression lookup
#-------------------------------------------------------------------------------------------------------

Execution: ./gene_expression_lookup_willerDiagram.sh -a /path/to/gene_annotations.txt -t /path/to/normalizedGTEx.tstat.txt -s /path/to/centroidSnps_final_grouped_snps_all.txt -n # -d # -e 0.# -i

Flags:
-a	: gene annotations file
-t	: file for normalized (rank-ordered) specific expression t-statistics
-n 	: number of genes to search per snp
-s 	: file containing snps to search
-d 	: window from snp within which to search for nearby genes (in kbp)
-e 	: threshold for high expression of a gene in a particular tissue
-i 	: include snp even if it has no nearby gene(s) within specified distance or if nearby gene(s) have missing expression t-statistics by finding the next nearest gene(s), even if these our outside the window

Output: 
	lookupResults_centroidSnps_final_grouped_snps_all.txt
	nearestGenes_centroidSnps_final_grouped_snps_all.txt

Notes: 
	Uses flags as inputs into expression_lookup_willerDiagram.py
	Must be in the same directory as expression_lookup_willerDiagram.py and nonOverlappingGenes.py
	Instead of the -i flag, can also use either -h or -z, which are described below: 
		-h: hide snp if it has no nearby gene(s) within specified distance or if nearby gene(s) have missing expression t-statistics by not including the snp in the expression output file
		-z: include snp even if it has no nearby gene within specified distance or if nearby gene(s) have missing expression t-statistics by having its expression vector be 0 (not highly expressed) for every tissue

#-------------------------------------------------------------------------------------------------------
# to create the ML table
#-------------------------------------------------------------------------------------------------------

Execution: python combineFeatures_newMLTables.py [1] [2]
	[1] = lookupResults_centroidSnps_final_grouped_snps_all.txt
	[2] = file of ML table with other desired features

Output:
	[2]_combined.txt 