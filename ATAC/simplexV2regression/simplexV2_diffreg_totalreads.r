library(randomForest)
library(caret)
library(doMC)
registerDoMC(cores = 10)
#library(glmnet)
# names(getModelInfo())

#data <- read.table("/scratch/sh8tv/Project/scATAC/Result/twoSideBias_correction/linear_regression_diffbias_reads/GM12878_CTCF/GM12878_CTCF_MC00128_GSM733752_combineALL.bed",header=T)
#paste0(outname,"_ATACsig_rf_CVpower.txt")
#outname <- "CTCF_MC00128_GSM733752"
#motif_cutoff <- 10000
#model_type <- "lm"

a<-commandArgs(T)
outname <- a[1]

data <- read.table(paste0("../combinesig_ALL/",outname,"_GC_peakov_combinesig.bed"))
colnames(data) <- c("chrm","start","end","motifseq","motifscore","strand","TFpeak","GC","ATACpeak","DNasepeak","TFsig","ATACsig","DNasesig")

motifscore <- log10(data[,"motifscore"])
TFpeak <- data[,"TFpeak"]
ATACpeak <- data[,"ATACpeak"]
DNasepeak <- data[,"DNasepeak"]
TFsig <- log10(data[,"TFsig"] + 1e-4)
ATACsig <- log10(data[,"ATACsig"] + 1e-4)
DNasesig <- log10(data[,"DNasesig"] + 1e-4)
GC <- data[,"GC"]


usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,ATACsig))

cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_ATACsig_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
rownames(powermat) <- c("rf")
write.table( powermat,file=paste0(outname,"_ATACsig_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)

cv_lm <- train(TFsig ~ .,data=usematrix,method="lm",trControl=trainControl(method="cv",number=5,verboseIter=T))
cv_step <- train(TFsig ~ .,data=usematrix,method="lmStepAIC",trControl=trainControl(method="cv",number=5,verboseIter=T))
cv_lasso <- train(TFsig ~ .,data=usematrix,method="lasso",trControl=trainControl(method="cv",number=5,verboseIter=T))

write.table( as.matrix(cv_lm$finalModel$coefficients),file=paste0(outname,"_ATACsig_rawlm_coeiff.txt"),sep="\t",quote=F,col.names=F,row.names=T)
write.table( as.matrix(cv_step$finalModel$coefficients),file=paste0(outname,"_ATACsig_step_coeiff.txt"),sep="\t",quote=F,col.names=F,row.names=T)

powermat <- rbind(apply(cv_lm$resample[,1:3],2,mean),
			apply(cv_step$resample[,1:3],2,mean),
			apply(cv_lasso$resample[,1:3],2,mean))
rownames(powermat) <- c("rawlm","step","lasso")
write.table( powermat,file=paste0(outname,"_ATACsig_lm_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)

save(cv_rf,file=paste0(outname,"_ATACsig_rfmodel.Rdata"))
save(cv_step,file=paste0(outname,"_ATACsig_stepmodel.Rdata"))


usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,DNasesig))

cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_DNasesig_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
rownames(powermat) <- c("rf")
write.table( powermat,file=paste0(outname,"_DNasesig_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)

cv_lm <- train(TFsig ~ .,data=usematrix,method="lm",trControl=trainControl(method="cv",number=5,verboseIter=T))
cv_step <- train(TFsig ~ .,data=usematrix,method="lmStepAIC",trControl=trainControl(method="cv",number=5,verboseIter=T))
cv_lasso <- train(TFsig ~ .,data=usematrix,method="lasso",trControl=trainControl(method="cv",number=5,verboseIter=T))

write.table( as.matrix(cv_lm$finalModel$coefficients),file=paste0(outname,"_DNasesig_rawlm_coeiff.txt"),sep="\t",quote=F,col.names=F,row.names=T)
write.table( as.matrix(cv_step$finalModel$coefficients),file=paste0(outname,"_DNasesig_step_coeiff.txt"),sep="\t",quote=F,col.names=F,row.names=T)

powermat <- rbind(apply(cv_lm$resample[,1:3],2,mean),
			apply(cv_step$resample[,1:3],2,mean),
			apply(cv_lasso$resample[,1:3],2,mean))
rownames(powermat) <- c("rawlm","step","lasso")
write.table( powermat,file=paste0(outname,"_DNasesig_lm_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)



save(cv_rf,file=paste0(outname,"_DNasesig_rfmodel.Rdata"))
save(cv_step,file=paste0(outname,"_DNasesig_stepmodel.Rdata"))



#usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,DNasesig))

#cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
#write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_DNasesig_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
#powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
#rownames(powermat) <- c("rf")
#write.table( powermat,file=paste0(outname,"_DNasesig_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)


#usematrix <- as.data.frame(cbind(TFsig,motifscore,GC,ATACsig,DNasesig))

#cv_rf <- train(TFsig ~ .,data=usematrix,method="rf",trControl=trainControl(method="cv",number=5,verboseIter=T))
#write.table( as.matrix(cv_rf$finalModel$importance),file=paste0(outname,"_bothsig_rf_importance.txt"),sep="\t",quote=F,col.names=F,row.names=T)
#powermat <- t(as.matrix(apply(cv_rf$resample[,1:3],2,mean)))
#rownames(powermat) <- c("rf")
#write.table( powermat,file=paste0(outname,"_bothsig_rf_CVpower.txt"),sep="\t",quote=F,col.names=T,row.names=T)
	


