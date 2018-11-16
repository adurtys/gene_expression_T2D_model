#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("3 arguments must be given (trait1, trait2, directory).n", call.=FALSE)
}

trait1 <- args[1]
trait2 <- args[2]
directory <- args[3]

#import data
trait1_train <-read.delim(file=paste(directory,"/",trait1,"_training_ML_table.txt",sep=""))
trait1_test  <-read.delim(file=paste(directory,"/",trait1,"_testing_ML_table.txt",sep=""))

trait2_train <-read.delim(file=paste(directory,"/",trait2,"_training_ML_table.txt",sep=""))
trait2_test  <-read.delim(file=paste(directory,"/",trait2,"_testing_ML_table.txt",sep=""))

#trait1_train <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer2\\T2D_like_training_ML_table.txt")
#trait1_test  <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer2\\T2D_like_testing_ML_table.txt")
#
#trait2_train <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer2\\lipid_training_ML_table.txt")
#trait2_test  <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer2\\lipid_testing_ML_table.txt")



#make datatables
trait1_train_data   <- trait1_train[,c(1,2)]
trait1_test_data    <- trait1_test[,c(1,2)]
trait2_train_data <- trait2_train[,c(1,2)]
trait2_test_data  <- trait2_test[,c(1,2)]

#trait1 model
trait1_train_data$trait1 <- (-2.41844289225346 + trait1_train$HepG2_ChIP.seq_POLR2AphosphoS5_ENCFF002CKW.bed*0.0360032900177372 + trait1_train$HepG2_ChIP.seq_MYBL2_ENCFF002CKR.bed*0.0013513548888414199 + trait1_train$E118.H3K27ac.gappedPeak*0.08704651878632902 + trait1_train$hepatocyte_H3K27Ac_DMSO.bed*0.265351986257085 + trait1_train$FAT.ADIP.NUC.EnhA*0.0102638702320791 + trait1_train$Hepg2.EnhA*0.0477263443760913 + trait1_train$Hepg2.Tx*0.0146133840269986 + trait1_train$hsap_HNF4A_hg19.bed*0.165938726233724 + trait1_train$GSE64233_MED1_V_final.bed*0.04809970248088599)
trait1_test_data$trait1 <- (-2.41844289225346 + trait1_test$HepG2_ChIP.seq_POLR2AphosphoS5_ENCFF002CKW.bed*0.0360032900177372 + trait1_test$HepG2_ChIP.seq_MYBL2_ENCFF002CKR.bed*0.0013513548888414199 + trait1_test$E118.H3K27ac.gappedPeak*0.08704651878632902 + trait1_test$hepatocyte_H3K27Ac_DMSO.bed*0.265351986257085 + trait1_test$FAT.ADIP.NUC.EnhA*0.0102638702320791 + trait1_test$Hepg2.EnhA*0.0477263443760913 + trait1_test$Hepg2.Tx*0.0146133840269986 + trait1_test$hsap_HNF4A_hg19.bed*0.165938726233724 + trait1_test$GSE64233_MED1_V_final.bed*0.04809970248088599)
trait2_train_data$trait1 <- (-2.41844289225346 + trait2_train$HepG2_ChIP.seq_POLR2AphosphoS5_ENCFF002CKW.bed*0.0360032900177372 + trait2_train$HepG2_ChIP.seq_MYBL2_ENCFF002CKR.bed*0.0013513548888414199 + trait2_train$E118.H3K27ac.gappedPeak*0.08704651878632902 + trait2_train$hepatocyte_H3K27Ac_DMSO.bed*0.265351986257085 + trait2_train$FAT.ADIP.NUC.EnhA*0.0102638702320791 + trait2_train$Hepg2.EnhA*0.0477263443760913 + trait2_train$Hepg2.Tx*0.0146133840269986 + trait2_train$hsap_HNF4A_hg19.bed*0.165938726233724 + trait2_train$GSE64233_MED1_V_final.bed*0.04809970248088599)
trait2_test_data$trait1 <- (-2.41844289225346 + trait2_test$HepG2_ChIP.seq_POLR2AphosphoS5_ENCFF002CKW.bed*0.0360032900177372 + trait2_test$HepG2_ChIP.seq_MYBL2_ENCFF002CKR.bed*0.0013513548888414199 + trait2_test$E118.H3K27ac.gappedPeak*0.08704651878632902 + trait2_test$hepatocyte_H3K27Ac_DMSO.bed*0.265351986257085 + trait2_test$FAT.ADIP.NUC.EnhA*0.0102638702320791 + trait2_test$Hepg2.EnhA*0.0477263443760913 + trait2_test$Hepg2.Tx*0.0146133840269986 + trait2_test$hsap_HNF4A_hg19.bed*0.165938726233724 + trait2_test$GSE64233_MED1_V_final.bed*0.04809970248088599)


#trait2 model
trait1_train_data$trait2 <- (-2.3528265626678397 + trait1_train$E087.H3K27ac.gappedPeak*0.017491268484867 + trait1_train$PE_Active_Enhancers_hg19.bed*0.0236999159838652 + trait1_train$Islets.stretchEnhancers.bed*0.115218392394746 + trait1_train$FOXA2*0.42479366225609794 + trait1_train$NKX2.2*0.19269221105605397 + trait1_train$PancIslt.EnhA*0.0589021324107632 + trait1_train$PDX1*0.34604711125045595)
trait1_test_data$trait2 <- (-2.3528265626678397 + trait1_test$E087.H3K27ac.gappedPeak*0.017491268484867 + trait1_test$PE_Active_Enhancers_hg19.bed*0.0236999159838652 + trait1_test$Islets.stretchEnhancers.bed*0.115218392394746 + trait1_test$FOXA2*0.42479366225609794 + trait1_test$NKX2.2*0.19269221105605397 + trait1_test$PancIslt.EnhA*0.0589021324107632 + trait1_test$PDX1*0.34604711125045595)
trait2_train_data$trait2 <- (-2.3528265626678397 + trait2_train$E087.H3K27ac.gappedPeak*0.017491268484867 + trait2_train$PE_Active_Enhancers_hg19.bed*0.0236999159838652 + trait2_train$Islets.stretchEnhancers.bed*0.115218392394746 + trait2_train$FOXA2*0.42479366225609794 + trait2_train$NKX2.2*0.19269221105605397 + trait2_train$PancIslt.EnhA*0.0589021324107632 + trait2_train$PDX1*0.34604711125045595)
trait2_test_data$trait2 <- (-2.3528265626678397 + trait2_test$E087.H3K27ac.gappedPeak*0.017491268484867 + trait2_test$PE_Active_Enhancers_hg19.bed*0.0236999159838652 + trait2_test$Islets.stretchEnhancers.bed*0.115218392394746 + trait2_test$FOXA2*0.42479366225609794 + trait2_test$NKX2.2*0.19269221105605397 + trait2_test$PancIslt.EnhA*0.0589021324107632 + trait2_test$PDX1*0.34604711125045595)

library(dplyr)
library(tidyr)

#loads ROCR package to calculate AUCs & ROC curves
library(ROCR)

####################################
#ROC plots!

#trait1_training
trait1_training.pred <- prediction(trait1_train_data$trait1,trait1_train_data$type)
trait1_training.perf  = performance(trait1_training.pred, measure = "tpr", x.measure = "fpr")
trait1_training.auc   = performance(trait1_training.pred, measure = "auc")

#trait2_training
trait2_training.pred <- prediction(trait2_train_data$trait2,trait2_train_data$type)
trait2_training.perf  = performance(trait2_training.pred, measure = "tpr", x.measure = "fpr")
trait2_training.auc   = performance(trait2_training.pred, measure = "auc")


#training AUCs
trait1_training.auc@y.values[[1]]
trait2_training.auc@y.values[[1]]

#training plots
pdf(paste(directory,"/",trait1,"_",trait2,"_training_ROC.pdf", sep = ""), width = 5, height = 5)

plot(trait2_training.perf, main=toupper("ROC for training"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(trait1_training.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait1,": ",round(trait1_training.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste(trait2,": ",round(trait2_training.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

dev.off()

#trait1 models on holdouts
trait1_model_trait1_testing.pred <-  prediction(trait1_test_data$trait1,trait1_test_data$type)
trait1_model_trait1_testing.perf  = performance(trait1_model_trait1_testing.pred, measure = "tpr", x.measure = "fpr")
trait1_model_trait1_testing.auc   = performance(trait1_model_trait1_testing.pred, measure = "auc")

trait1_model_trait2_testing.pred <-  prediction(trait2_test_data$trait1,trait2_test_data$type)
trait1_model_trait2_testing.perf  = performance(trait1_model_trait2_testing.pred, measure = "tpr", x.measure = "fpr")
trait1_model_trait2_testing.auc   = performance(trait1_model_trait2_testing.pred, measure = "auc")


trait1_model_trait1_testing.auc@y.values
trait1_model_trait2_testing.auc@y.values

#trait2 models on holdouts
trait2_model_trait1_testing.pred <-  prediction(trait1_test_data$trait2,trait1_test_data$type)
trait2_model_trait1_testing.perf  = performance(trait2_model_trait1_testing.pred, measure = "tpr", x.measure = "fpr")
trait2_model_trait1_testing.auc   = performance(trait2_model_trait1_testing.pred, measure = "auc")

trait2_model_trait2_testing.pred <-  prediction(trait2_test_data$trait2,trait2_test_data$type)
trait2_model_trait2_testing.perf  = performance(trait2_model_trait2_testing.pred, measure = "tpr", x.measure = "fpr")
trait2_model_trait2_testing.auc   = performance(trait2_model_trait2_testing.pred, measure = "auc")


trait2_model_trait1_testing.auc@y.values
trait2_model_trait2_testing.auc@y.values



#plot it!
pdf(paste(directory,"/",trait1,"_",trait2,"_testing_ROC.pdf", sep = ""), width = 10, height = 5)

par(mfrow=c(1,2))
plot(trait1_model_trait2_testing.perf, main=toupper(paste(trait1,"model on holdouts")), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(trait1_model_trait1_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait1,":",round(trait1_model_trait1_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste(trait2,":",round(trait1_model_trait2_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

plot(trait2_model_trait2_testing.perf, main=toupper(paste(trait2,"model on holdouts")), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(trait2_model_trait1_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait1,":",round(trait2_model_trait1_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste(trait2,":",round(trait2_model_trait2_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

dev.off()

#compare to CADD, GWAVA, DeepSEA
trait1_GCD <- read.delim(file=paste(directory,"/",trait1,"_testing_scores_table.txt",sep=""))

trait2_GCD <- read.delim(file=paste(directory,"/",trait2,"_testing_scores_table.txt",sep=""))

#remember to invert DeepSEA data
trait1_GCD$deepsea_flip <- (1-trait1_GCD$DeepSEA)
trait2_GCD$deepsea_flip <- (1-trait2_GCD$DeepSEA)



#ROC plots!

#CADD on holdouts
CADD_trait1_testing.pred <-  prediction(trait1_GCD$CADD,trait1_GCD$type)
CADD_trait1_testing.perf  = performance(CADD_trait1_testing.pred, measure = "tpr", x.measure = "fpr")
CADD_trait1_testing.auc   = performance(CADD_trait1_testing.pred, measure = "auc")

CADD_trait2_testing.pred <-  prediction(trait2_GCD$CADD,trait2_GCD$type)
CADD_trait2_testing.perf  = performance(CADD_trait2_testing.pred, measure = "tpr", x.measure = "fpr")
CADD_trait2_testing.auc   = performance(CADD_trait2_testing.pred, measure = "auc")

CADD_trait1_testing.auc@y.values
CADD_trait2_testing.auc@y.values

#GWAVA on holdouts
GWAVA_trait1_testing.pred <-  prediction(trait1_GCD$GWAVA,trait1_GCD$type)
GWAVA_trait1_testing.perf  = performance(GWAVA_trait1_testing.pred, measure = "tpr", x.measure = "fpr")
GWAVA_trait1_testing.auc   = performance(GWAVA_trait1_testing.pred, measure = "auc")

GWAVA_trait2_testing.pred <-  prediction(trait2_GCD$GWAVA,trait2_GCD$type)
GWAVA_trait2_testing.perf  = performance(GWAVA_trait2_testing.pred, measure = "tpr", x.measure = "fpr")
GWAVA_trait2_testing.auc   = performance(GWAVA_trait2_testing.pred, measure = "auc")

GWAVA_trait1_testing.auc@y.values
GWAVA_trait2_testing.auc@y.values

#DeepSEA on holdouts
DeepSEA_trait1_testing.pred <-  prediction(trait1_GCD$deepsea_flip,trait1_GCD$type)
DeepSEA_trait1_testing.perf  = performance(DeepSEA_trait1_testing.pred, measure = "tpr", x.measure = "fpr")
DeepSEA_trait1_testing.auc   = performance(DeepSEA_trait1_testing.pred, measure = "auc")

DeepSEA_trait2_testing.pred <-  prediction(trait2_GCD$deepsea_flip,trait2_GCD$type)
DeepSEA_trait2_testing.perf  = performance(DeepSEA_trait2_testing.pred, measure = "tpr", x.measure = "fpr")
DeepSEA_trait2_testing.auc   = performance(DeepSEA_trait2_testing.pred, measure = "auc")

DeepSEA_trait1_testing.auc@y.values
DeepSEA_trait2_testing.auc@y.values


#plot it!
pdf(paste(directory,"/testing_ROC_wCGD.pdf", sep = ""), width = 10, height = 5)

par(mfrow=c(1,2))
plot(trait2_model_trait1_testing.perf, main=toupper(paste("models on",trait1,"holdouts")), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(trait1_model_trait1_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
plot(DeepSEA_trait1_testing.perf, col = "#e6ab02", lwd = 2)
par(new=TRUE)
plot(CADD_trait1_testing.perf, col = "#66a61e", lwd = 2)
par(new=TRUE)
plot(GWAVA_trait1_testing.perf, col = "#e7298a", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait1,":",round(trait1_model_trait1_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste(trait2,":",round(trait2_model_trait1_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)
text(1, .2, paste("CADD:",round(CADD_trait1_testing.auc@y.values[[1]],2),sep=""), col = "#66a61e", adj = 1)
text(1, .1, paste("GWAVA:",round(GWAVA_trait1_testing.auc@y.values[[1]],2),sep=""), col = "#e7298a", adj = 1)
text(1, 0, paste("DeepSEA:",round(DeepSEA_trait1_testing.auc@y.values[[1]],2),sep=""), col = "#e6ab02", adj = 1)

plot(trait2_model_trait2_testing.perf, main=toupper(paste("models on",trait2,"holdouts")), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(trait1_model_trait2_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
plot(DeepSEA_trait2_testing.perf, col = "#e6ab02", lwd = 2)
par(new=TRUE)
plot(CADD_trait2_testing.perf, col = "#66a61e", lwd = 2)
par(new=TRUE)
plot(GWAVA_trait2_testing.perf, col = "#e7298a", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait1,":",round(trait1_model_trait2_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste(trait2,":",round(trait2_model_trait2_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)
text(1, .2, paste("CADD:",round(CADD_trait2_testing.auc@y.values[[1]],2),sep=""), col = "#66a61e", adj = 1)
text(1, .1, paste("GWAVA:",round(GWAVA_trait2_testing.auc@y.values[[1]],2),sep=""), col = "#e7298a", adj = 1)
text(1, 0, paste("DeepSEA:",round(DeepSEA_trait2_testing.auc@y.values[[1]],2),sep=""), col = "#e6ab02", adj = 1)

dev.off()


#and for training data:
trait1_tr_GCD <- read.delim(file=paste(directory,"/",trait1,"_training_scores_table.txt",sep=""))

trait2_tr_GCD <- read.delim(file=paste(directory,"/",trait2,"_training_scores_table.txt",sep=""))

#invert DeepSea
trait1_tr_GCD$deepsea_flip <- (1-trait1_tr_GCD$DeepSEA)
trait2_tr_GCD$deepsea_flip <- (1-trait2_tr_GCD$DeepSEA)


#CADD
CADD_trait1_training.pred <-  prediction(trait1_tr_GCD$CADD,trait1_tr_GCD$type)
CADD_trait1_training.perf  = performance(CADD_trait1_training.pred, measure = "tpr", x.measure = "fpr")
CADD_trait1_training.auc   = performance(CADD_trait1_training.pred, measure = "auc")

CADD_trait2_training.pred <-  prediction(trait2_tr_GCD$CADD,trait2_tr_GCD$type)
CADD_trait2_training.perf  = performance(CADD_trait2_training.pred, measure = "tpr", x.measure = "fpr")
CADD_trait2_training.auc   = performance(CADD_trait2_training.pred, measure = "auc")

CADD_trait1_training.auc@y.values
CADD_trait2_training.auc@y.values

#GWAVA
GWAVA_trait1_training.pred <-  prediction(trait1_tr_GCD$GWAVA,trait1_tr_GCD$type)
GWAVA_trait1_training.perf  = performance(GWAVA_trait1_training.pred, measure = "tpr", x.measure = "fpr")
GWAVA_trait1_training.auc   = performance(GWAVA_trait1_training.pred, measure = "auc")

GWAVA_trait2_training.pred <-  prediction(trait2_tr_GCD$GWAVA,trait2_tr_GCD$type)
GWAVA_trait2_training.perf  = performance(GWAVA_trait2_training.pred, measure = "tpr", x.measure = "fpr")
GWAVA_trait2_training.auc   = performance(GWAVA_trait2_training.pred, measure = "auc")

GWAVA_trait1_training.auc@y.values
GWAVA_trait2_training.auc@y.values

#DeepSEA
DeepSEA_trait1_training.pred <-  prediction(trait1_tr_GCD$deepsea_flip,trait1_tr_GCD$type)
DeepSEA_trait1_training.perf  = performance(DeepSEA_trait1_training.pred, measure = "tpr", x.measure = "fpr")
DeepSEA_trait1_training.auc   = performance(DeepSEA_trait1_training.pred, measure = "auc")

DeepSEA_trait2_training.pred <-  prediction(trait2_tr_GCD$deepsea_flip,trait2_tr_GCD$type)
DeepSEA_trait2_training.perf  = performance(DeepSEA_trait2_training.pred, measure = "tpr", x.measure = "fpr")
DeepSEA_trait2_training.auc   = performance(DeepSEA_trait2_training.pred, measure = "auc")

DeepSEA_trait1_training.auc@y.values
DeepSEA_trait2_training.auc@y.values


pdf(paste(directory,"/training_ROC_wCGD.pdf", sep = ""), width = 10, height = 5)

par(mfrow=c(1,2))
plot(trait1_training.perf, main=toupper(paste(trait1,"training")), col = "forestgreen", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(DeepSEA_trait1_training.perf, col = "#e6ab02", lwd = 2)
par(new=TRUE)
plot(CADD_trait1_training.perf, col = "#66a61e", lwd = 2)
par(new=TRUE)
plot(GWAVA_trait1_training.perf, col = "#e7298a", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait1,":",round(trait1_training.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste("CADD:",round(CADD_trait1_training.auc@y.values[[1]],2),sep=""), col = "#66a61e", adj = 1)
text(1, .2, paste("GWAVA:",round(GWAVA_trait1_training.auc@y.values[[1]],2),sep=""), col = "#e7298a", adj = 1)
text(1, .1, paste("DeepSEA:",round(DeepSEA_trait1_training.auc@y.values[[1]],2),sep=""), col = "#e6ab02", adj = 1)

plot(trait2_training.perf, main=toupper(paste(trait2,"training")), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(DeepSEA_trait2_training.perf, col = "#e6ab02", lwd = 2)
par(new=TRUE)
plot(CADD_trait2_training.perf, col = "#66a61e", lwd = 2)
par(new=TRUE)
plot(GWAVA_trait2_training.perf, col = "#e7298a", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
text(1, .5, "AUC", adj = 1)
text(1, .4, paste(trait2,":",round(trait2_training.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)
text(1, .3, paste("CADD:",round(CADD_trait2_training.auc@y.values[[1]],2),sep=""), col = "#66a61e", adj = 1)
text(1, .2, paste("GWAVA:",round(GWAVA_trait2_training.auc@y.values[[1]],2),sep=""), col = "#e7298a", adj = 1)
text(1, .1, paste("DeepSEA:",round(DeepSEA_trait2_training.auc@y.values[[1]],2),sep=""), col = "#e6ab02", adj = 1)

dev.off()
