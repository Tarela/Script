library(randomForest)
library(caret)
#library(doParallel)
#registerDoParallel(cores=10)
#library(doMC)
#registerDoMC(cores = 10)
#library(glmnet)
# names(getModelInfo())

#data <- read.table("/scratch/sh8tv/Project/scATAC/Result/twoSideBias_correction/linear_regression_diffbias_reads/GM12878_CTCF/GM12878_CTCF_MC00128_GSM733752_combineALL.bed",header=T)
#paste0(outname,"_",motif_cutoff,"_plus_rf_CVpower.txt")
#outname <- "CTCF_MC00128_GSM733752"
#motif_cutoff <- 10000
#model_type <- "lm"

a<-commandArgs(T)
outname <- a[1]

data <- read.table(paste0(outname,"_simplexV2output.bed"))
colnames(data) <- c("chrm","start","end","motifseq","motifscore","strand","TFpeak","ATACpeak","TFsig","ATACsig","GC",
					"plus.intercept", paste0("plus.mono.",seq(1:24)), paste0("plus.di.",seq(1:63)),
					"minus.intercept", paste0("minus.mono.",seq(1:24)), paste0("minus.di.",seq(1:63))
					)

motifscore <- log10(data[,"motifscore"])
TFpeak <- data[,"TFpeak"]
ATACpeak <- data[,"ATACpeak"]
TFsig <- log10(data[,"TFsig"] + 1e-4)
ATACsig <- log10(data[,"ATACsig"] + 1e-4)
GC <- data[,"GC"]

plus_param <- data[,12:99]
minus_param <- data[,100:187]
cb_param <- plus_param + minus_param
colnames(cb_param) <- c("cb.intercept", paste0("cb.mono.",seq(1:24)), paste0("cb.di.",seq(1:63)))



usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,plus_param))

cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_plus_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
#rownames(powermat) <- c("rf")
write.table( powermat,file=paste0(outname,"_plus_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=F)
save(cv_rf,file=paste0(outname,"_plus_rfmodel.Rdata"))

cv_lm <- train(TFsig ~ .,data=usematrix,method="lm",trControl=trainControl(method="cv",number=5,verboseIter=T))
cv_step <- train(TFsig ~ .,data=usematrix,method="lmStepAIC",trControl=trainControl(method="cv",number=5,verboseIter=T))
cv_lasso <- train(TFsig ~ .,data=usematrix,method="lasso",trControl=trainControl(method="cv",number=5,verboseIter=T))

write.table( as.matrix(cv_lm$finalModel$coefficients),file=paste0(outname,"_plus_rawlm_coeiff.txt"),sep="\t",quote=F,col.names=F,row.names=T)
write.table( as.matrix(cv_step$finalModel$coefficients),file=paste0(outname,"_plus_step_coeiff.txt"),sep="\t",quote=F,col.names=F,row.names=T)

powermat <- rbind(apply(cv_lm$resample[,1:3],2,mean),
			apply(cv_step$resample[,1:3],2,mean),
			apply(cv_lasso$resample[,1:3],2,mean))
rownames(powermat) <- c("rawlm","step","lasso")
write.table( powermat,file=paste0(outname,"_plus_lm_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)
save(cv_step,file=paste0(outname,"_plus_stepmodel.Rdata"))




#usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,minus_param))

#cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
#write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_","_minus_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
#powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
#rownames(powermat) <- c("rf")
#write.table( powermat,file=paste0(outname,"_","_minus_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)


#usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,cb_param))

#cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
#write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_","_combine_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
#powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
#rownames(powermat) <- c("rf")
#write.table( powermat,file=paste0(outname,"_","_combine_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)
	


