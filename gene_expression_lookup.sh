# Date Created: 9 July 2018
# Date Last Modified: 1 August 2018
# Execution: ./gene_expression_lookup_final.sh (TODO WITH FLAGSS!)
# Description: Shell script for tissue expression lookup of genes near specified snps. First, ensures that both the 
# 	GENCODE GTF file and GTEx t-statistics file have been processed. Then analyzes expression of genes closest to snps
# 	listed in the snp file.
# 	Note: files used in the flags must be contained in the working directory. 

#!/usr/bin/env bash

onlyProtCoding=false
processMissingSnps="I"

while getopts :a:d:e:ghin:prs:t:z option
do
	case "$option" in
		a)
			# GENCODE GTF file has already been processed
			geneAnnotationsFile="$OPTARG"
			echo "Gene annotations filename: $geneAnnotationsFile"
			;;
		d)
			distanceFromSnp="$OPTARG"
			echo "Distance from snp to search for genes: $distanceFromSnp"
			;;
		e)
			expressionThreshold="$OPTARG"
			echo "Threshold for high expression: $expressionThreshold"
			;;
		g)
			# process GENCODE GTF file to create a gene annotations file
			echo "Process GENCODE GTF file."
			processGencode=true
			;;
		h)
			# hide snps whose nearest gene doesn't have tissue expression t-stats
			echo "If nearest gene to snp doesn't have tissue expression t-stats, exclude snp from expression vector."
			processMissingSnps="H"
			;;
		i)
			# include snps whose nearest gene doesn't have tissue expression t-stats by using t-stat for next nearest gene
			echo "If nearest gene to snp doesn't have tissue expression t-stats, use t-stats from next nearest gene that has t-stats."
			# processMissingSnps="I"
			;;
		n)
			# search for this number of nearest genes
			numNearestGenesToSearch="$OPTARG"
			echo "Number of nearest genes to search from snp: $numNearestGenesToSearch"
			;;
		p)
			# only include protein-coding genes in gene_annotations file
			echo "Only include protein-coding genes in gene annotations file."
			onlyProtCoding=true
			;;

		r)
			# rank-order t-stats for tissue expression
			echo "Normalize t-statistics in GTEx file."
			python ./tstat_normalization.py
			tStatFile="./normalizedGTEx.tstat.txt"
			;;
		s)
			snpFile="$OPTARG"
			echo "File containing snps: $snpFile"

			# create output filenames
			outFile=lookupResults_$snpFile
			nearestGenesFile=nearestGenes_$snpFile
			;;
		t)
			# GTEx t-statistics file has been normalized
			tStatFile="$OPTARG"
			echo "Normalized t-stat filename: $tStatFile"
			;;
		z)
			# include snps whose nearest gene doesn't have tissue expression t-stats by using 0 for all tissues in the expression vector
			echo "If nearest gene to snp doesn't have tissue expression t-stats, use 0 in tissue expression vector."
			processMissingSnps="Z"
			;;
		:)
			echo "ERROR: Option -$OPTARG requires an argument."
			exit 1
			;;
		\?) 
			echo "ERROR: Invalid option: -$OPTARG"
			exit 1
			;;
	esac
done

# determine whether / how to process GENCODE file
if [ "$processGencode" = true ]
then
	echo "Have to process GENCODE file."
	if [ "$onlyProtCoding" = true ]
	then
		echo "Processing GENCODE file. Only including protein-coding genes."
		python ./GTF_processing.py gencode.v19.annotation.gtf genesInGTEx_v6p.txt $tStatFile true
	else [ "$onlyProtCoding" = false ]
		echo "Processing GENCODE file. Including all genes in tstat file in addition to protein-coding genes."
		python ./GTF_processing.py gencode.v19.annotation.gtf genesInGTEx_v6p.txt $tStatFile false
	fi

	geneAnnotationsFile="./gene_annotations.txt"
else
	echo "GENCODE file has already been processed."
fi

# TODO: figure out whether or not to normalize GTEx tstat file

# conduct expression lookup
python ./expression_lookup.py $snpFile $geneAnnotationsFile $tStatFile $numNearestGenesToSearch $distanceFromSnp $expressionThreshold $outFile $nearestGenesFile $processMissingSnps

