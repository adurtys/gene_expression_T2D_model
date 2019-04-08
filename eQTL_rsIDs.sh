# Date Created: 8 April 2019
# Date Last Modified: 8 April 2019
# Execution: ./eQTL_rsIDs.sh
# Description: obtain the rsID corresponding to each GWAS snp in LD with a variant

#!/usr/bin/env bash
for var in ./isLD.eQTL_variants_reformatted.txt
do
	grep $var ./all_gregor_snplist_rsIDs.txt
done
