library(ROCR)
#### new50 100
AR50new<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/newform/AR_FP_50_100bp_peak.bed")
AR50newData <- AR50new[which(AR50new[,7]<0),]
#### new50-100
CTCF50new<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/newform/CTCF_FP_50_100bp_peak.bed")
CTCF50newData <- CTCF50new[which(CTCF50new[,7]<0),]
#### new25U50U
AR25U50Unew<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/newform/AR_FP_25U50U_peak.bed")
AR25U50UnewData <- AR25U50Unew[which(AR25U50Unew[,7]<0),]

#### new25U50U
CTCF25U50Unew<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/newform/CTCF_FP_25U50U_peak.bed")
CTCF50newData <- CTCF25U50Unew[which(CTCF25U50Unew[,7]<0),]
#### newSEPE
ARSEPEnew<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/newform/AR_FP_SEPE_peak.bed")
ARSEPEnewData <- ARSEPEnew[which(ARSEPEnew[,7]<0),]

#### newSEPE
CTCFSEPEnew<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/newform/CTCF_FP_SEPE_peak.bed")
CTCFSEPEnewData <- CTCFSEPEnew[which(CTCFSEPEnew[,7]<0),]

#### old 50-100
AR50old<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/oldform/AR_FP_old_peak.bed")
AR50oldData<-AR50oldData[which(AR50oldData[,7]!="no"),]

CTCF50old<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/oldform/CTCF_FP_peak.bed")
CTCF50oldData<-CTCF50old

#### oldSEPE
ARSEPEold<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/oldform/AR_FP_SEPE_peak.bed")
ARSEPEoldData<-ARSEPEold

CTCFSEPEold<-read.table("~/Desktop/DNAse/scan_footprint/binding_FPscore/oldform/CTCF_FP_SEPE_peak.bed")
CTCFSEPEoldData<-CTCFSEPEold


##ROC curve
#AR 50-100
par(mfrow=c(2,2),mar=c(4,4,2,2))
pred<-prediction(as.numeric(as.vector(AR50oldData[,7])),as.numeric(as.vector(AR50oldData[,10])))
perf<-performance(pred,"tpr","fpr")
perfA<-performance(pred,"auc")
PA<-perfA@"y.values"[[1]]
plot(perf,col="darkblue",main="AR 50-100bp")
par(new=T)
pred<-prediction(as.numeric(as.vector(AR50newData[,7])),as.numeric(as.vector(AR50newData[,10])))
perfnew<-performance(pred,"tpr","fpr")
perfB<-performance(pred,"auc")
PB<-perfB@"y.values"[[1]]
plot(perfnew,col="darkred",axes=F,xlab="",ylab="",main="")
legend("topleft",c(paste0("Old formula AUC=",round(PA,3)),paste0("New formula AUC=",round(PB,3))),col=c("darkblue","darkred"),lwd=3)
#CTCF 50-100
pred<-prediction(as.numeric(as.vector(CTCF50oldData[,7])),as.numeric(as.vector(CTCF50oldData[,10])))
perf<-performance(pred,"tpr","fpr")
perfA<-performance(pred,"auc")
PA<-perfA@"y.values"[[1]]
plot(perf,col="darkblue",main="CTCF 50-100bp")
par(new=T)
pred<-prediction(as.numeric(as.vector(CTCF50newData[,7])),as.numeric(as.vector(CTCF50newData[,10])))
perfnew<-performance(pred,"tpr","fpr")
perfB<-performance(pred,"auc")
PB<-perfB@"y.values"[[1]]
plot(perfnew,col="darkred",axes=F,xlab="",ylab="",main="")
legend("topleft",c(paste0("Old formula AUC=",round(PA,3)),paste0("New formula AUC=",round(PB,3))),col=c("darkblue","darkred"),lwd=3)
#AR SEPE
pred<-prediction(as.numeric(as.vector(ARSEPEoldData[,7])),as.numeric(as.vector(ARSEPEoldData[,10])))
perf<-performance(pred,"tpr","fpr")
perfA<-performance(pred,"auc")
PA<-perfA@"y.values"[[1]]
plot(perf,col="darkblue",main="AR SEPE")
par(new=T)
pred<-prediction(as.numeric(as.vector(ARSEPEnewData[,7])),as.numeric(as.vector(ARSEPEnewData[,10])))
perfnew<-performance(pred,"tpr","fpr")
perfB<-performance(pred,"auc")
PB<-perfB@"y.values"[[1]]
plot(perfnew,col="darkred",axes=F,xlab="",ylab="",main="")
legend("topleft",c(paste0("Old formula AUC=",round(PA,3)),paste0("New formula AUC=",round(PB,3))),col=c("darkblue","darkred"),lwd=3)
#CTCF SEPE
pred<-prediction(as.numeric(as.vector(CTCFSEPEoldData[,7])),as.numeric(as.vector(CTCFSEPEoldData[,10])))
perf<-performance(pred,"tpr","fpr")
perfA<-performance(pred,"auc")
PA<-perfA@"y.values"[[1]]
plot(perf,col="darkblue",main="CTCF SEPEbp")
par(new=T)
pred<-prediction(as.numeric(as.vector(CTCFSEPEnewData[,7])),as.numeric(as.vector(CTCFSEPEnewData[,10])))
perfnew<-performance(pred,"tpr","fpr")
perfB<-performance(pred,"auc")
PB<-perfB@"y.values"[[1]]
plot(perfnew,col="darkred",axes=F,xlab="",ylab="",main="")
legend("topleft",c(paste0("Old formula AUC=",round(PA,3)),paste0("New formula AUC=",round(PB,3))),col=c("darkblue","darkred"),lwd=3)



####Density plot
par(mfrow=c(4,2),mar=c(4,4,2,2))
### 50-100bp old fomula
#AR
plot(density(as.numeric(as.vector(AR50oldData[which(AR50oldData[,10]==0),7])),bw=0.1),xlim=c(-5,5),ylim=c(0,3),col="darkred",xlab="FPscore",ylab="Density",main="Old AR 25~50U 50~100bp DHT+Veh")
par(new=T)
plot(density(as.numeric(as.vector(AR50oldData[which(AR50oldData[,10]!=0),7])),bw=0.1),xlim=c(-5,5),ylim=c(0,3),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)
#CTCF
plot(density(as.numeric(as.vector(CTCF50oldData[which(CTCF50oldData[,10]==0),7]))),xlim=c(-50,50),ylim=c(0,0.1),col="darkred",xlab="FPscore",ylab="Density",main="Old CTCF 25~50U 50~100bp DHT+Veh")
par(new=T)
plot(density(as.numeric(as.vector(CTCF50oldData[which(CTCF50oldData[,10]!=0),7]))),xlim=c(-50,50),ylim=c(0,0.1),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)
### 50-100bp new fomula
#AR
plot(density(as.numeric(as.vector(AR50newData[which(AR50newData[,10]==0),7])),bw=0.05),xlim=c(-1.5,0),ylim=c(0,3),col="darkred",xlab="FPscore",ylab="Density",main="AR 25~50U 50~100bp DHT+Veh new fomula ")
par(new=T)
plot(density(as.numeric(as.vector(AR50newData[which(AR50newData[,10]!=0),7])),bw=0.05),xlim=c(-1.5,0),ylim=c(0,3.5),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)
#CTCF
plot(density(as.numeric(as.vector(CTCF50newData[which(CTCF50newData[,10]==0),7])),bw=0.05),xlim=c(-2.5,0),ylim=c(0,3),col="darkred",xlab="FPscore",ylab="Density",main="New CTCF 25~50U 50~100bp DHT+Veh")
par(new=T)
plot(density(as.numeric(as.vector(CTCF50newData[which(CTCF50newData[,10]!=0),7])),bw=0.05),xlim=c(-2.5,0),ylim=c(0,3),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)

#par(mfrow=c(2,2),mar=c(4,4,2,2))
### 50-100bp old fomula
#AR
plot(density(as.numeric(as.vector(ARSEPEoldData[which(ARSEPEoldData[,10]==0),7])),bw=0.1),xlim=c(-5,5),ylim=c(0,1),col="darkred",xlab="FPscore",ylab="Density",main="Old AR all SEPE")
par(new=T)
plot(density(as.numeric(as.vector(ARSEPEoldData[which(ARSEPEoldData[,10]!=0),7])),bw=0.1),xlim=c(-5,5),ylim=c(0,1),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)
#CTCF
plot(density(as.numeric(as.vector(CTCFSEPEoldData[which(CTCFSEPEoldData[,10]==0),7]))),xlim=c(-50,50),ylim=c(0,0.1),col="darkred",xlab="FPscore",ylab="Density",main="Old CTCF all SEPE")
par(new=T)
plot(density(as.numeric(as.vector(CTCFSEPEoldData[which(CTCFSEPEoldData[,10]!=0),7]))),xlim=c(-50,50),ylim=c(0,0.1),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)
### 50-100bp new fomula
#AR
plot(density(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]==0),7]))[which(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]==0),7]))<0)],bw=0.05),xlim=c(-1.5,0),ylim=c(0,3),col="darkred",xlab="FPscore",ylab="Density",main="AR all SEPE new fomula")
par(new=T)
plot(density(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]!=0),7]))[which(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]!=0),7]))<0)],bw=0.05),xlim=c(-1.5,0),ylim=c(0,3.5),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)
#CTCF
plot(density(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]==0),7]))[which(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]==0),7]))<0)],bw=0.05),xlim=c(-2.5,0),ylim=c(0,3),col="darkred",xlab="FPscore",ylab="Density",main="New CTCF all SEPE")
par(new=T)
plot(density(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]!=0),7]))[which(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]!=0),7]))<0)],bw=0.05),xlim=c(-2.5,0),ylim=c(0,3),col="darkblue",axes=F,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)

#### method compare
png(filename="method_compare_peak_nopeak.png",width=900,height=900)
par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]==0),7])),as.numeric(as.vector(ARSEPEoldData[which(ARSEPEoldData[,10]==0),7])),xlim=c(-1.5,0),ylim=c(-7,7),xlab="SEPEnew",ylab="SEPE old",main="AR , fomula result compare on SEPE",col="#FF0000",pch=".")
par(new=T)
plot(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]!=0),7])),as.numeric(as.vector(ARSEPEoldData[which(ARSEPEoldData[,10]!=0),7])),col="#0000FF",pch=".",xlim=c(-1.5,0),ylim=c(-7,7),axes=T,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)


plot(as.numeric(as.vector(AR50newData[which(AR50newData[,10]==0),7])),as.numeric(as.vector(AR50oldData[which(AR50oldData[,10]==0),7])),xlim=c(-1.5,0),ylim=c(-7,7),xlab="50new",ylab="50 old",main="AR , fomula result compare on 50-100bp",col="#FF0000",pch=".")
par(new=T)
plot(as.numeric(as.vector(AR50newData[which(AR50newData[,10]!=0),7])),as.numeric(as.vector(AR50oldData[which(AR50oldData[,10]!=0),7])),col="#0000FF",pch=".",xlim=c(-1.5,0),ylim=c(-7,7),axes=T,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)

plot(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]==0),7])),as.numeric(as.vector(CTCFSEPEoldData[which(CTCFSEPEoldData[,10]==0),7])),xlim=c(-1.5,0),ylim=c(-70,50),xlab="SEPEnew",ylab="SEPE old",main="CTCF , fomula result compare on SEPE",col="#FF000050",pch=".")
par(new=T)
plot(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]!=0),7])),as.numeric(as.vector(CTCFSEPEoldData[which(CTCFSEPEoldData[,10]!=0),7])),col="#0000FF50",pch=".",xlim=c(-1.5,0),ylim=c(-70,50),axes=T,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)


plot(as.numeric(as.vector(CTCF50newData[which(CTCF50newData[,10]==0),7])),CTCF50old_nopeak,xlim=c(-1.5,0),ylim=c(-70,50),xlab="50new",ylab="50 old",main="CTCF , fomula result compare on 50-100bp",col="#FF000050",pch=".")
par(new=T)
plot(as.numeric(as.vector(CTCF50newData[which(CTCF50newData[,10]!=0),7])),CTCF50old_peak,col="#0000FF50",pch=".",xlim=c(-1.5,0),ylim=c(-70,50),axes=T,xlab="",ylab="",main="")
legend("topleft",c("overlapChippeak","no Chippeak"),col=c("darkblue","darkred"),lwd=3)

dev.off()

#### method compare
png(filename="method_compare_sep.png",width=1600,height=900)
par(mfcol=c(2,4),mar=c(4,4,2,2))
smoothScatter(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]==0),7])),as.numeric(as.vector(ARSEPEoldData[which(ARSEPEoldData[,10]==0),7])),main="AR SEPE noPeak",ylim=c(-7,7),xlim=c(-1,0))
smoothScatter(as.numeric(as.vector(ARSEPEnewData[which(ARSEPEnewData[,10]!=0),7])),as.numeric(as.vector(ARSEPEoldData[which(ARSEPEoldData[,10]!=0),7])),main="AR SEPE Peak",ylim=c(-7,7),xlim=c(-1,0))

smoothScatter(as.numeric(as.vector(AR50newData[which(AR50newData[,10]==0),7])),as.numeric(as.vector(AR50oldData[which(AR50oldData[,10]==0),7])),main="AR 50-100bp noPeak",ylim=c(-7,7),xlim=c(-1,0))
smoothScatter(as.numeric(as.vector(AR50newData[which(AR50newData[,10]!=0),7])),as.numeric(as.vector(AR50oldData[which(AR50oldData[,10]!=0),7])),main="AR 50-100bp Peak",ylim=c(-7,7),xlim=c(-1,0))

smoothScatter(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]==0),7])),as.numeric(as.vector(CTCFSEPEoldData[which(CTCFSEPEoldData[,10]==0),7])),main="CTCF SEPE noPeak",ylim=c(-70,50),xlim=c(-1.5,0))
smoothScatter(as.numeric(as.vector(CTCFSEPEnewData[which(CTCFSEPEnewData[,10]!=0),7])),as.numeric(as.vector(CTCFSEPEoldData[which(CTCFSEPEoldData[,10]!=0),7])),main="CTCF SEPE Peak",ylim=c(-70,50),xlim=c(-1.5,0))

smoothScatter(as.numeric(as.vector(CTCF50newData[which(CTCF50newData[,10]==0),7])),CTCF50old_nopeak,main="CTCF 50-100bp noPeak",ylim=c(-70,50),xlim=c(-1.5,0))
smoothScatter(as.numeric(as.vector(CTCF50newData[which(CTCF50newData[,10]!=0),7])),CTCF50old_peak,main="CTCF 50-100bp Peak",ylim=c(-70,50),xlim=c(-1.5,0))
dev.off()


### logistic regression
x1<-as.numeric(as.vector(CTCF50newData[,7]))###FPscore
x2<-sqrt(as.numeric(as.vector(CTCF50newData[,8])))###DNase count
y<-as.numeric(as.vector(CTCF50newData[,10]))###label

all.data <- data.frame(cbind(y, x1, x2))
train.data.index <- sample(1:nrow(all.data), size = floor(nrow(all.data) * 0.8))
train.data <- all.data[train.data.index,]
###
test.data <- all.data[-(train.data.index),]
###
logistic.model_12 <- glm(y ~ x1 + x2, family = binomial, data = train.data)
logistic.model_1 <- glm(y ~ x1, family = binomial, data = train.data)
logistic.model_2 <- glm(y ~ x2, family = binomial, data = train.data)

logistic.pred_12 <- predict(logistic.model_12, test.data[,2:3]) #test.data
logistic.pred_1 <- logistic.model_1$coeff[1]+logistic.model_1$coeff[2]*test.data[,2] 
logistic.pred_2 <- logistic.model_2$coeff[1]+logistic.model_2$coeff[2]*test.data[,3] 
###logistic.model$coeff
###logistic.model$coeff[1] + logistic.model$coeff[2] * test.data$x1 + logistic.model$coeff[3] * test.data$x2
logistic.p_12 <- inv.logit(logistic.pred_12)
logistic.p_1 <- inv.logit(logistic.pred_1)
logistic.p_2 <- inv.logit(logistic.pred_2)

comp_12 <- prediction(logistic.p_12, test.data[,1])
perf_12 <- performance(comp_12, "tpr", "fpr")
perfauc_12 <- performance(comp_12, "auc")
comp_1 <- prediction(logistic.p_1, test.data[,1])
perf_1 <- performance(comp_1, "tpr", "fpr")
perfauc_1 <- performance(comp_1, "auc")
comp_2 <- prediction(logistic.p_2, test.data[,1])
perf_2 <- performance(comp_2, "tpr", "fpr")
perfauc_2 <- performance(comp_2, "auc")

plot(perf_12,col="darkred",lwd=2,xlab="False positive rate",ylab="Ture positive rate",main="ROC curve, single time cross validation")
par(new=T)
plot(perf_1,col="darkblue",lwd=2,xlab="",ylab="",main="",axes=F)
par(new=T)
plot(perf_2,col="darkgreen",lwd=2,xlab="",ylab="",main="",axes=F)
legend("bottomright",c("Tagcount+FPscore","Tagcount","FPscore"),col=c("darkred","darkgreen","darkblue"),lwd=4)
plot(perf)
plot(firth.perf)

###1000次， a,b1,b2取平均值
auc12<-c()
auc1<-c()
auc2<-c()
for(i in 1:1000){
train.data.index <- sample(1:nrow(all.data), size = floor(nrow(all.data) * 0.8))
train.data <- all.data[train.data.index,]
###
test.data <- all.data[-(train.data.index),]
###
logistic.model_12 <- glm(y ~ x1 + x2, family = binomial, data = train.data)
logistic.model_1 <- glm(y ~ x1, family = binomial, data = train.data)
logistic.model_2 <- glm(y ~ x2, family = binomial, data = train.data)

logistic.pred_12 <- predict(logistic.model_12, test.data[,2:3]) #test.data
logistic.pred_1 <- logistic.model_1$coeff[1]+logistic.model_1$coeff[2]*test.data[,2] 
logistic.pred_2 <- logistic.model_2$coeff[1]+logistic.model_2$coeff[2]*test.data[,3] 
###logistic.model$coeff
###logistic.model$coeff[1] + logistic.model$coeff[2] * test.data$x1 + logistic.model$coeff[3] * test.data$x2
logistic.p_12 <- inv.logit(logistic.pred_12)
logistic.p_1 <- inv.logit(logistic.pred_1)
logistic.p_2 <- inv.logit(logistic.pred_2)

comp_12 <- prediction(logistic.p_12, test.data[,1])
perfauc_12 <- performance(comp_12, "auc")@"y.values"[[1]]
comp_1 <- prediction(logistic.p_1, test.data[,1])
perfauc_1 <- performance(comp_1, "auc")@"y.values"[[1]]
comp_2 <- prediction(logistic.p_2, test.data[,1])
perfauc_2 <- performance(comp_2, "auc")@"y.values"[[1]]
auc12<-c(auc12,perfauc_12)
auc1<-c(auc1,perfauc_1)
auc2<-c(auc2,perfauc_2)
}
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(density(auc12),col="darkred",xlab="Area under curve(AUC) of ROC",ylab="Density",main="Distribution of AUC in 240 times cross validation",lwd=2,ylim=c(0,360),xlim=c(0.7,1))
par(new=T)
plot(density(auc1),col="darkblue",xlab="",ylab="",main="",axes=F,lwd=2,ylim=c(0,360),xlim=c(0.7,1))
par(new=T)
plot(density(auc2),col="darkgreen",xlab="",ylab="",main="",axes=F,lwd=2,ylim=c(0,360),xlim=c(0.7,1))
legend("topleft",c("Tagcount+FPscore","Tagcount","FPscore"),col=c("darkred","darkgreen","darkblue"),lwd=4)

boxplot(auc12,auc1,auc2,col=c("red","blue","green"),main="Boxplot of AUC",ylab="AUC")
legend("topleft",c("Tagcount+FPscore","Tagcount","FPscore"),col=c("darkred","darkgreen","darkblue"),lwd=4)

###library(brglm)
library(brglm)
firth.model <- brglm(y ~ x1 + x2, methods = "brglm", data = train.data)
firth.pred <- predict(firth.model, test.data[,2:3])
firth.p <- inv.logit(firth.pred)
firth.comp <- prediction(firth.p, test.data[,1])
firth.perf <- performance(firth.comp, "tpr", "fpr")
