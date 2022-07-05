library(ROCR)
a<-commandArgs(T)
infile <- a[1]
outfile <- a[2]
data<-read.table(infile)
label <- data[,7]
label[which(label > 1)] <- 1
TC <- data[,8]
FOS <- data[,9]
FPmin <- data[,10]
FPmax <- data[,11]
FPmax[FPmax == 0] <- 1
FPmax<-  -log(FPmax)


pred_TC<-prediction(as.numeric(as.vector(TC)),as.numeric(as.vector(label)))
perf_TC<-performance(pred_TC,"tpr","fpr")
perfA_TC<-performance(pred_TC,"auc",fpr.stop=0.05)

pred_FOS<-prediction(as.numeric(as.vector(FOS)),as.numeric(as.vector(label)))
perf_FOS<-performance(pred_FOS,"tpr","fpr")
perfA_FOS<-performance(pred_FOS,"auc",fpr.stop=0.05)

pred_FP<-prediction(as.numeric(as.vector(FPmax)),as.numeric(as.vector(label)))
perf_FP<-performance(pred_FP,"tpr","fpr")
perfA_FP<-performance(pred_FP,"auc",fpr.stop=0.05)

## combine TC and FP using logistic regression
y<-as.numeric(as.vector(label))
x1<-as.numeric(as.vector(TC))
x2<-as.numeric(as.vector(FPmax))

all.data <- data.frame(cbind(y, x1, x2))
train.data.index <- sample(1:nrow(all.data), size = floor(nrow(all.data) * 0.8))
train.data <- all.data[train.data.index,]
test.data <- all.data[-(train.data.index),]

logistic.model_12 <- glm(y ~ x1 + x2, family = binomial, data = train.data)
logistic.pred_12 <- predict(logistic.model_12, test.data[,2:3])
logistic.p_12 <- inv.logit(logistic.pred_12)

pred_12 <- prediction(logistic.p_12, test.data[,1])
perf_12 <- performance(pred_12, "tpr", "fpr")
perfauc_12 <- performance(pred_12, "auc",fpr.stop=0.05)


pdf(paste0(outfile,".pdf"))
plot(perf_TC,col="red")
par(new=T)
plot(perf_FOS,xlab="",ylab="",axes=F,col="blue")
par(new=T)
plot(perf_FP,xlab="",ylab="",axes=F,col="green")
par(new=T)
plot(perf_12,xlab="",ylab="",axes=F,col="purple")

legend("topleft",c(paste0("TC  AUC0.05 = ",round(perfA_TC@y.values[[1]],4)),paste0("ownFOS  AUC0.05 = ",round(perfA_FOS@y.values[[1]],4)),paste0("FPregion  AUC0.05 = ",round(perfA_FP@y.values[[1]],4)),paste0("TC+FP  AUC0.05 = ",round(perfauc_12@y.values[[1]],4))),col=c("red","blue","green","purple"),lwd=4)
dev.off()

print(c(a[2],round(perfA_TC@y.values[[1]],4),round(perfA_FOS@y.values[[1]],4),round(perfA_FP@y.values[[1]],4),round(perfauc_12@y.values[[1]],4)))


