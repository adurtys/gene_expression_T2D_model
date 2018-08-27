# Date Created: 3 August 2018
# Date Last Modified: 20 August 2018
# Execution: ./centroid_snp_lookup.sh
# Description: TODO

#!/usr/bin/env bash

for i in final_grouped_snps_chr*.txt
do
	echo "Finding centroid snps for groups in $i"
	python centroidPositions.py $i
done

# create one file with centroid snps for all chromosomes
echo "Combining individual chromosome files into one file, 'centroidSnps.txt'"
for j in centroidSnps_final_grouped_snps_chr*.txt
do
	if [ $j = 0 ]
	then
		echo "Creating centroidSnps.txt"
		cat $j > centroidSnps.txt
	else
		cat $j >> centroidSnps.txt
	fi
done
