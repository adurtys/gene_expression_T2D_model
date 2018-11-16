#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("four arguments must be given (MLtable file, random number file, analysis name, out directory).n", call.=FALSE)
  }
#read in data
MLtable<-read.table(file=args[1], header= TRUE)
random_numbers <- read.table(file=args[2], header = FALSE)

#make datatable
datatable <- MLtable[,c(1,2)]
  
#slice and dice matrices
MLtable <- data.matrix(MLtable)
y <- MLtable[,"type"]
x <- MLtable[,-1:-2]



##keeping this for future but do not have this data right now
#loads ROCR package and creates predictions for various scores
#
#CADD.pred <- prediction(datatable$CADD_score,datatable$snp_group)
#CADD.perf = performance(CADD.pred, measure = "tpr", x.measure = "fpr")
#CADD.auc = performance(CADD.pred, measure = "auc")
#
#GWAVA_TSS.pred <- prediction(datatable$GWAVA_tss_score,datatable$snp_group)
#GWAVA_TSS.perf = performance(GWAVA_TSS.pred, measure = "tpr", x.measure = "fpr")
#GWAVA_TSS.auc = performance(GWAVA_TSS.pred, measure = "auc")
#
#FATHMM_MKL_NC.pred <- prediction(datatable$FATHMM.MKL_NC,datatable$snp_group)
#FATHMM_MKL_NC.perf = performance(FATHMM_MKL_NC.pred, measure = "tpr", x.measure = "fpr")
#FATHMM_MKL_NC.auc = performance(FATHMM_MKL_NC.pred, measure = "auc")
#
#PhastCon_100.pred <- prediction(datatable$PhastCon_100way,datatable$snp_group)
#PhastCon_100.perf = performance(PhastCon_100.pred, measure = "tpr", x.measure = "fpr")
#PhastCon_100.auc = performance(PhastCon_100.pred, measure = "auc")
#
##print established loci prediction AUCs
#CADD.auc@y.values
#GWAVA_TSS.auc@y.values
#FATHMM_MKL_NC.auc@y.values
#PhastCon_100.auc@y.values
#

#loads glmnet & ROCR packages and creates 15 10fold CV models
library(glmnet)
library(ROCR)

count = 0
for (i in random_numbers$V1) {  
  count = count+1
  #set folds
  set.seed(i)
  foldid=sample(1:10,size=length(y),replace=TRUE)
  
  #build model
  model_out = cv.glmnet(x, y, foldid = foldid, family = "binomial", alpha = 1, type.measure="auc")
  coef_out <- coef(model_out, s = "lambda.1se")
  predict_out <- predict(model_out, newx = x, s = "lambda.1se")
  
  #add to datatable
  columnname <- paste ("model_",count,sep="")
  datatable[[columnname]] <- predict_out
  
  #use ROCR to calc AUC
  model.pred <- prediction(datatable[[columnname]],datatable$type)
  model.perf = performance(model.pred, measure = "tpr", x.measure = "fpr")
  model.auc = performance(model.pred, measure = "auc")
  AUC = model.auc@y.values
  print(paste("model_",count," AUC=",AUC,sep="") )
  #save coef to dataframe
  coef_mat <- as(coef_out, "matrix")
  coef_mat <- rbind(coef_mat,"AUC"=AUC[[1]][1])
  #if first time through, save as dataframe
  if (count == 1) {
    coef_df <- as.data.frame(coef_mat)
    colnames(coef_df) <- columnname
    
  }else { #for rest just append to dataframe
    coef_df[[columnname]] <- coef_mat
    
  }
  
}
#print coefficients to file
write.table(coef_df, paste(args[4],args[3],"_coefficients.txt",sep=''), quote = FALSE, sep="\t", row.names = TRUE)

#print datatable to file
write.table(datatable, paste(args[4],args[3],"_models.txt",sep=''), quote = FALSE, sep="\t", row.names = FALSE)
