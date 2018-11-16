# gene_expression_T2D_model
ENTIRE PIPELINE
---------------
 
1. create gene annotations file for autosomal genes (from GENCODE) contained in GTEx --> currently working with all genes, regardless of whether or not they are protein coding (TODO: analysis with only protein-coding genes?)
 
python GTF_processing.py gencode.v19.annotation.gtf GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct GTEx.tstat.tsv false
 
2. normalize t-statistics for tissue expression (rank-orders them)
 
python tstat_normalization.py GTEx.tstat.tsv
 
3. create centroid snp list for all grouped snps on all chromosomes
 
./centroid_snp_lookup.sh
 
4. (MANUAL) create .txt files with only snpGroup and snpType (index or control) for each of the four categories (lipid training, lipid testing, T2D-like training, T2D-like testing)
 
cut -f1,2 category_ML_table_filename.txt> [category]Snps.txt
* should have 4 files total - one for each category
 
5. create file containing the categories for each grouped snp
 
python snpsInModel.py lipidTestingSnps.txt lipidTrainingSnps.txt T2DLikeTestingSnps.txt T2DLikeTrainingSnps.txt
 
 
(only have to do the above once - same outputs can be used in all model runs)
-----------------------------------------------------------------------------
 
6. conduct expression lookup
 
./gene_expression_lookup.sh [OPTIONS]
 
in our case (if above three steps were already conducted):
./gene_expression_lookup.sh -a gene_annotations.txt -t normalizedGTEx.tstat.txt -s centroidSnps.txt -n # -d 1000 -e 0.# -i
NOTE: searches for nearest genes within 1 mbp (because -d 1000)
 
7. create the ML tables containing only gene expression features
 
python expression_ML_table.py snpTypes.txt lookupResults_centroidSnps.txt
* should have 4 files total - one for each category (but only have to run once)
 
8. create ML tables that have all features (combines expression ML tables to ML tables lacking expression)
 
python combineFeatures.py snpCategory_expression_ML_table.txt snpCategory_noExpression_ML_table.txt
* should run four times - once for each category - to get 4 combined ML table output files, total
 
9. (MANUAL) generate 15 random seeds for lipid model and 15 random seeds for T2D-like model (6-7 digit number)
 
(below steps involve scripts written by Kim)
--------------------------------------------
 
10. run LASSO to build models based off of training set
 
Rscript --vanilla path/to/10foldCV_LASSO_glmnet.R path/to/training/ML_table.txt path/to/random_seeds.txt prefix path/to/output/directory(./ for current wd)
* each test should run six times - 3  expression model types (noExpression, onlyExpression, and combined) for 2 traits (lipid and T2D)
 
11. parse coefficients, removing rows that are all 0
 
perl /project/voight_ML/lorenzk/V2/scripts/Parse_coef_output.pl snpCategory_expressionModel_coefficients.txt snpCategory_expressionModel_coefficients_parsed.txt
* run six times total - once for each model type
 
12. (MANUAL) scan parsed coefficients file to select features - make trait_expressionModel_model.txt (first row has intercept, then two columns - first is feature name, second is coefficient, where features are listed by descending order of coefficients)
 
* done six times total - once for each model type
 
13. convert trait_expressionModel_model.txt to proper format for R script
 
perl /project/voight_ML/lorenzk/V2/scripts/print_R_model.pl path/to/T2D_training_ML_table.txt,path/to/T2D_testing_ML_table.txt,path/to/lipid_training_ML_table.txt,path/to/lipid_testing_ML_table.txt path/to/trait_expressionModel_model.txt > trait_expressionModel_model_Rformat.txt
* done six times total - once for each model type
 
14. (MANUAL) paste these R-format vectors into the local R script (testing_sets_cross_pred_combined.R) and create ROC plots in RStudio
 
source("testing_sets_cross_pred_expressionModel.R", echo = TRUE)
 
---------------------------------------------------------------------------------------------------------------------------------------
SCRIPTS (that I wrote)
GTF_processing.py --> gene_annotations.txt
tstat_normalization.py --> normalizedGTEx.tstat.txt
centroidPositions.py && centroid_snp_lookup.sh --> centroidSnps.txt
snpsInModel.py --> snpTypes.txt
nonOverlappingGenes.py --> genesWithoutTstat.txt
expression_lookup.py && gene_expression_lookup.sh ---> nearestGenes_centroidSnps.txt && lookupResults_centroidSnps.txt
expression_ML_table.py --> snpCategory_expression_ML_table.txt
combineFeatures.py --> snpCategory_combined_ML_table.txt
------------------------------------------------------------------------------------------------------------------------------------------

GTF_processing.py: creates a text file containing annotations (start and end locations) for autosomal genes

argv1: GENCODE GTF filename
argv2: GTEx v6p filename
argv3: filename containing tissue expression t-statistics
argv4: boolean for whether gene annotations file should include only protein-coding genes ("true"), or if it should include all genes that have tissue expression t-statistics, regardless of whether or not they are protein-coding ("false")
Expected Output: gene_annotations.txt

column[0]: ENSGID (no decimal places following the ID)
column[1]: gene name
column[2]: chromosome number (format: "chr:#")
column[3]: gene start location
column[4]: gene end location
column[5]: feature type (whether the gene is "protein_coding" or not)
python GTF_processing.py gencode.v19.annotation.gtf GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct GTEx.tstat.tsv false

wc -l gene_annotations.txt --> 25214 genes
python GTF_processing.py gencode.v19.annotation.gtf GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct GTEx.tstat.tsv true

wc -l gene_annotations.txt --> 18296 genes
Notes:

File excludes non-autosomal genes (all genes on chromosomes X, Y, and M)
File excludes genes not in GTEx v6p (because these do not have tissue expression t-statistics)
372 genes in the annotations file do not have tissue expression t-statistics (these are saved locally in "genesInGTEx_v6p_noTstat.txt")
tstat_normalization.py: rank-orders the tissue expression t-statistics file

argv1: filename containing tissue expression t-statistics
Expected Output: normalizedGTEx.tstat.txt

organized the same way as original t-statistics file (header line, first column contains ENSGID)
python tstat_normalization.py GTEx.tstat.tsv

wc -l normalizedGTEx.tstat.txt --> 24843 lines (incl. header)
centroidPositions.py: determines centroid snp for a list of all snps in each group of snps associated by LD on a given chromosome

argv1: groupedSnpsFilename
column[0]: group number
column[1]: snp - in format chr#:snpLocation
Expected Output: centroidSnps_groupedSnpsFilename.txt

Notes:

centroid snp location was determined by taking the average of the minimum and maximum snp locations for the entire group of snps
centroid_snp_lookup.sh: finds centroid snps for each of the chromosomes in final_grouped_snps directory

Expected Output:

centroidSnps.txt
centroidSnps_groupedSnpsFilename.txt for each of the chromosomes in the directory
./centroid_snp_lookup.sh

wc -l centroidSnps.txt --> 10237 snps
Notes: should be run inside the final_grouped_snps directory, which should contain grouped snp files for all (22 autosomal) chromosomes
snpsInModel.py: creates file containing characterization of each snp group with respect to type (index or control) and category (lipid testing, lipid training, T2D-like testing, T2D-like training)

argv1: filename for file containing lipid testing snps and their type
argv2: filename for file containing lipid training snps and their type
argv3: filename for file containing T2D-like testing snps and their type
argv4: filename for file containing T2D-like training snps and their type
Expected Output: snpTypes.txt

column[0]: snpGroup
column[1]: snpType --> index or control
column[2]: snpCategory --> lipidTesting, lipidTraining, T2DLikeTesting, or T2DLikeTraining
Notes:

Four input files were created manually, using cut command on Kim's ML tables
Some snps (~20) are in multiple categories
python snpsInModel.py lipidTestingSnps.txt lipidTrainingSnps.txt T2DLikeTestingSnps.txt T2DLikeTrainingSnps.txt

 wc -l snpTypes.txt --> 8097 (including header line, so 8096 snp groups in total)
nonOverlappingGenes.py: finds genes contained in gene annotations file but that do not have tissue expression t-statistics

argv1: filename for file containing tissue expression t-statistics
argv2: filename for gene annotations file
Expected Output: genesWithoutTstat.txt

python nonOverlappingGenes.py GTEx.tstat.txt gene_annotations.txt

wc -l genesWithoutTstat.txt --> 372 genes in gene annotations file that do not have tissue expression t-statistics
Notes: 
Test on when gene_annotations.txt includes all genes with tissue-expression t-statistics (not only the protein-coding ones)
expression_lookup.py: conducts the expression lookup for the nearest gene to each snp in the input file

argv1: filename for file containing snps to search
argv2: gene annotations filename
argv3: filename for file containing rank-orders of the tissue expression t-statistics
argv4: number of nearest genes to the snp for which to conduct the tissue expression lookup
argv5: distance from the snp for which to look for nearby genes
argv6: the percent of top ranks of t-statistics that should be considered "highly expressed" for each tissue
argv7: filename for output file containing tissue expression vectors for nearest numGenes to each snp
argv8: filename for output file containing nearest numGenes to each snp
argv9: flag for how to process snps whose nearest gene has missing tissue-expression data, either because there is no gene within the specified distance from the snp or because any or all of the nearest numGenes don't have tissue expression t-statistics
Expected Output:

nearestGenes[File].txt
column[0]: snp grouping
column[1]: snp
column[2]: nearest gene(s) depending on numGenes and whether snp is equidistant from genes
outputFile.txt
gene_expression_lookup.sh: shell script to run lookup pipeline

Possible Flags:
-a [geneAnnotationsFilename]: gene annotations file has already been created (filename used as argv2 in expression_lookup.py)
-g: create gene annotations file --> runs GTF_processing.py with arguments hardcoded (and in same directory)
-t [normalizedTstatFilename]: t-statistics file has already been normalized (filename used argv3 in expression_lookup.py)
-r: normalize (rank-order) the t-statistics file --> runs t-stat_normalization.py
-n [numGenes]: number of genes to search per snp (used as argv4 in expression_lookup.py)
-s [snpFilename]: file containing snps to search (used as argv1 in expression_lookup.py)
-d [distanceFromSnp]: distance from snp to include in search for nearby genes (used as argv5 in expression_loookup.py)
-e [expressionThreshold]: threshold for "high expression" of a gene in a particular tissue (used as argv6 in expression_lookup.py)
-i: include snp with no nearby gene (within specified distance) or with missing tissue expression t-statistics by finding the next nearest gene that meets this criteria
-h: hide snp with no nearby gene (within specified distance) or with missing tissue expression t-statistics by not including the snp in the tissue expression output file
-z: include snp with no nearby gene (within specified distance) or with missing tissue expression t-statistics by having its expression vector be 0 (not highly expressed) for every tissue
Expected Output:

nearestGenes_snpFilename
lookupResults_snpFilename
expression_ML_table.py: creates ML tables containing only gene expression features

argv1: filename for file containing snp type for each group of snps (index or control)
argv2: filename for file containing grouped snp gene expression lookup results (determined by the centroid snp's expression for each group)
Expected Output: 4 files - one for each category (lipid training, lipid testing, T2D-like training, T2D-like testing)

snpCategory_expression_ML_table.txt
python expression_ML_table.py snpTypes.txt lookupResults_centroidSnps.txt

Notes:

categories are hard-coded right now (TODO: make this more general?)
combineFeatures.py: combines ML tables that lack expression features and contain only expression features, respectively

argv1: filename for file containing ML table containing only gene expression features
argv2: filename for file containing ML table that lacks gene expression features
Expected Output: snpCategory_combined_ML_table.txt
