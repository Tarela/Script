library(ROCR)
a<-commandArgs(T)
infile <- a[1]
outfile <- a[2]
data<-read.table(infile)
label <- data[,7]
label[which(label > 1)] <- 1
TC <- data[,8]
FOS <- data[,9]

pred_TC<-prediction(as.numeric(as.vector(TC)),as.numeric(as.vector(label)))
perf_TC<-performance(pred_TC,"tpr","fpr")
perfA_TC<-performance(pred_TC,"auc")

pred_FOS<-prediction(as.numeric(as.vector(FOS)),as.numeric(as.vector(label)))
perf_FOS<-performance(pred_FOS,"tpr","fpr")
perfA_FOS<-performance(pred_FOS,"auc")

pdf(paste0(outfile,".pdf"))
plot(perf_TC,col="red")
par(new=T)
plot(perf_FOS,xlab="",ylab="",axes=F,col="blue")
legend("topleft",c(paste0("TC  AUC = ",round(perfA_TC@y.values[[1]],4)),paste0("FOS  AUC = ",round(perfA_FOS@y.values[[1]],4))),col=c("red","blue"),lwd=4)
dev.off()

print(c(a[2],round(perfA_TC@y.values[[1]],4),round(perfA_FOS@y.values[[1]],4)))


