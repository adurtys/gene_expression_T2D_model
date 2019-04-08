# Date Created: 8 April 2019
# Date Last Modified: 8 April 2019
# Execution: python isLD_reformatted.py [1]
# argv[1] = isLD.eQTL_variants.txt
# Description: reformats the eQTL variants in order to determine rsIDs when creating isQTL matrix
# Run Time: ~1 sec

#!/usr/bin/env python
import sys

# read in command-line arguments
isLD_variant_id_filename = sys.argv[1]

isLD_variant_id_file = open(isLD_variant_id_filename, 'r')

variant_dict = {} # key = variant_id; value = reformatted_variant_id (chr#:pos_first_ref_base)
for line in isLD_variant_id_file:
	variant_id = line.rstrip('\r\n')
	reformatted_variant_id = variant_id
	if variant_id not in variant_dict:
		variant_info = variant_id.split("_")
		reformatted_variant_id = "chr" + variant_info[0] + ":" + variant_info[1]

		variant_dict[variant_id] = reformatted_variant_id

isLD_variant_id_file.close()

tab = "\t"
newline = "\n"

output = ""

for variant in variant_dict:
	output += variant + tab + variant_dict[variant] + newline

outfilename = "isLD.eQTL_variants_reformatted.txt"
outfile = open(outfilename, 'w')
outfile.write(output)
outfile.close()