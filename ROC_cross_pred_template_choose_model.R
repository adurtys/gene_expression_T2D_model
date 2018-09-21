#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("4 arguments must be given (trait1, trait2, model, directory).n", call.=FALSE)
}

model <- args[3]
trait1 <- args[1]
trait2 <- args[2]
directory <- args[4]

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

###insert models here

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
pdf(paste(directory,"/",trait1,"_",trait2,"_",model,"_training_ROC.pdf", sep = ""), width = 5, height = 5)

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
pdf(paste(directory,"/",trait1,"_",trait2,"_",model,"_testing_ROC.pdf", sep = ""), width = 10, height = 5)

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
