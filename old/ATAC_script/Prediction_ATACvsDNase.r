library(ROCR)
a<-commandArgs(T)
infile <- a[1]
outfile <- a[2]
data<-read.table(infile)
label <- data[,7]
label[which(label > 1)] <- 1
ATAC100 <- data[,8]
ATAC247 <- data[,9]
ATACall <- data[,10]
DNase <- data[,11]
pred_ATAC100<-prediction(as.numeric(as.vector(ATAC100)),as.numeric(as.vector(label)))
perf_ATAC100<-performance(pred_ATAC100,"tpr","fpr")
perfA_ATAC100<-performance(pred_ATAC100,"auc")

pred_ATAC247<-prediction(as.numeric(as.vector(ATAC247)),as.numeric(as.vector(label)))
perf_ATAC247<-performance(pred_ATAC247,"tpr","fpr")
perfA_ATAC247<-performance(pred_ATAC247,"auc")

pred_ATACall<-prediction(as.numeric(as.vector(ATACall)),as.numeric(as.vector(label)))
perf_ATACall<-performance(pred_ATACall,"tpr","fpr")
perfA_ATACall<-performance(pred_ATACall,"auc")

pred_DNase<-prediction(as.numeric(as.vector(DNase)),as.numeric(as.vector(label)))
perf_DNase<-performance(pred_DNase,"tpr","fpr")
perfA_DNase<-performance(pred_DNase,"auc")
rain <- rainbow(4)
pdf(file=paste0(outfile,".pdf"))
plot(perf_ATAC100,col=rain[1])
par(new=T)
plot(perf_ATAC247,col=rain[2],axes=F,xlab="",ylab="")
par(new=T)
plot(perf_ATACall,col=rain[3],xlab="",ylab="",axes=F)
par(new=T)
plot(perf_DNase,col=rain[4],xlab="",ylab="",axes=F)

legend("bottomright",c(paste0("ATAC100 AUC = ",round(perfA_ATAC100@y.values[[1]],4)),paste0("ATAC247 AUC = ",round(perfA_ATAC247@y.values[[1]],4)),paste0("ATACall AUC = ",round(perfA_ATACall@y.values[[1]],4)),paste0("DNase AUC = ",round(perfA_DNase@y.values[[1]],4))),col=rain,lwd=4,bty="n")
dev.off()

print(c(a[2],round(perfA_ATAC100@y.values[[1]],4),round(perfA_ATAC247@y.values[[1]],4),round(perfA_ATACall@y.values[[1]],4),round(perfA_DNase@y.values[[1]],4)))

