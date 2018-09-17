import sys
import numpy as np
import pandas
from StringIO import StringIO
import os
import pdb
import gzip
import argparse
import subprocess
import os.path
import string
import re
from time import sleep
from collections import defaultdict

###################################################
##program goals
## take in information (model type, model names, R template, & directory) that describe output of 10foldCV_LASSO_glmnet.R script
## this version is set up to process models for 2 traits at once
## parse coefficient information (remove coefficients with 0)
## return summary of model information to terminal & ask user to select model to use
## insert model coefficient information into R script to generate ROC plots for training & testing data
## run R script to obtain ROC plots
###################################################

parser = argparse.ArgumentParser()
parser.add_argument("-n", help="model names",type=str)
parser.add_argument("-t", help="R script template",type=str)
parser.add_argument("-o", help="directory to put stuff into",type=str)
parser.add_argument("-m", help="model type",type=str)


args = parser.parse_args()

Name = args.n
namelist = Name.split(",")
template = args.t
Out = args.o
modeltype = args.m

###################################################
##definitions
###################################################

#send to command line
def runstuff(cmd):
	subprocess.call(cmd, shell=True)
	#print("submitted "+cmd)
	
#send to command line, return output as list of lines
def get_out(cmd):
	stuff = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#print("returned "+cmd)
	stdout = []
	while True:
		l = stuff.stdout.readline()
		line = l.rstrip()
		if line == '' and stuff.poll() != None:
			break
		else:
			stdout.append(line)
	return stdout

###################################################
##select 1 of 15 models from glmnet output
###################################################
	
#step through each model name separately
modellist = []
for i in namelist:
	###################################################
	##run coef parse script & consolidate coef information
	###################################################
	Infile = Out+"/"+modeltype+"_"+i+"_coefficients.txt"
	parsed_coef = Out+"/"+modeltype+"_"+i+"_coefficients_parsed.txt"
	runstuff("perl /project/voight_ML/lorenzk/V2/scripts/Parse_coef_output.pl "+Infile+" "+parsed_coef)
	
	#read back in parsed coef table
	coef_table = pandas.read_table(parsed_coef, index_col='features')
	
	#get info to summarize for each of 15 models
	models = list(coef_table)
	infohash = defaultdict(list)
	AUClist=[]
	for x in models:
		xAUC = round(coef_table.loc['AUC',x],3)
		xfeat = str(coef_table.loc['#_selected_features',x])
		xint = str(round(coef_table.loc['(Intercept)',x],3))
		
		AUClist.append(xAUC)
		
		uniq_model = str(xAUC)+"\t"+xfeat+"\t"+xint
		infohash[uniq_model].append(x)

	AUC_med = np.median(AUClist)
	AUC_mean = round(np.mean(AUClist),3)
	#return model summaries to user:
	print("\n"+modeltype+" "+i+" models summary data")
	print("AUC median: "+str(AUC_med)+"\nAUC mean: "+str(AUC_mean))
	print("\nAUC\t#feat\tintercept\t#models\tmodels")
	for key, value in infohash.iteritems():
		count = len(value)
		print (key+"\t\t"+str(count)+"\t"+" ".join(str(x) for x in value))
	print(" ")
	###################################################
	##prompt for model selection & then pull that column
	##print model to file & run R parse script
	###################################################
	#this will prompt, comment out to streamline testing
	model = raw_input('please select model:')
	
	#model = "model_1"
	
	#add chosen model to list
	modellist.append(model)
	
	model_coef_filepath = Out+"/"+modeltype+"_"+i+"_"+model+".txt"
	
	#removing coef with 0 value
	select = coef_table[coef_table[model] != 0]
	
	#take just column with model we want
	select2 = select[model].copy()
	
	#remove #features row and AUC row
	select3 = select2.drop(['AUC', '#_selected_features'])
	
	#print to file
	select3.to_csv(model_coef_filepath, sep="\t")
	
	#run R formatting script with this model
	R_formatted_model = Out+"/"+modeltype+"_"+i+"_"+model+"_R.txt"
	runstuff("perl /project/voight_ML/lorenzk/V2/scripts/print_R_model.pl trait1_train,trait1_test,trait2_train,trait2_test "+model_coef_filepath+" > "+R_formatted_model)
	
###################################################
##read in R template & insert models generated above into script
##save as named R script & run
###################################################
outR_path = Out+"/"+modeltype+"_"+namelist[0]+"_"+modellist[0]+"_"+namelist[1]+"_"+modellist[1]+"_cross_pred_ROC.R"
outR = open(outR_path, 'w')

for l in open(template,'r'):
	line = l.rstrip()
	#if insert models line, open files and insert
		
	if "###insert models here" in line:
		#start with first trait
		R_formatted_model = Out+"/"+modeltype+"_"+namelist[0]+"_"+modellist[0]+"_R.txt"
		outR.write("#trait1 model\n")
		for q in open(R_formatted_model, 'r'):
			qline = q.rstrip()
			
			#pull out trait# from line
			matching = re.search('\+ (.+?)\$',qline).group(1)
			
			outR.write(matching+"_data$trait1 <- "+qline+"\n")
		
		
		#and then do the same for trait2 (notice the print statement is different!)
		R_formatted_model = Out+"/"+modeltype+"_"+namelist[1]+"_"+modellist[1]+"_R.txt"
		outR.write("\n\n#trait2 model\n")
		for q in open(R_formatted_model, 'r'):
			qline = q.rstrip()
			
			matching = re.search('\+ (.+?)\$',qline).group(1)
			
			outR.write(matching+"_data$trait2 <- "+qline+"\n")
		
	else:
		outR.write(line+"\n")
		

#run the R script!
runstuff("bsub -o "+Out+"/ROC_Rscript.out -q voight_normal Rscript --vanilla "+outR_path+" "+namelist[0]+" "+namelist[1]+" "+modeltype+" "+Out)