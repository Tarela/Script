a<-commandArgs(T)

library(hydroGOF)     
data<-read.table(a[1])
bias_span <-  a[3]
tag_span <- a[4]
fname <- a[2]

outrmse<-c()
outslope<-c()
outcc<-c()

png(paste0(fname,"_mean.png"),width=960,height=960)
par(mar=c(6,6,6,2))
plot(log10(data[,1]),log10(data[,2]),cex.axis=2.5,cex.lab=2.5,cex.main=3,pch=".",xlab=paste0("log10 observed average cut within +-",bias_span,"bp"),ylab=paste0("log10 predicted cut with +-",bias_span,"bp average bias"),main=paste0("bias window +-",bias_span,"bp\ntotal reads interval +-",tag_span,"bp"))

linR <- lm(log10(data[,2]+1e-5)~log10(data[,1]+1e-5))
abline(linR,lwd=3,col="red")

legend("bottomright",
        c(  paste0("rmse(linear) = ",round(rmse(data[,1]+1e-5,data[,2]+1e-5),4)),
            paste0("slope(log) = ",round(as.numeric(linR[[1]][2]),4)),
            paste0("C.C(log) = ",round(cor(log10(data[,1]+1e-5),log10(data[,2]+1e-5)),4)) ),
        bty="n", cex=3  )
dev.off()
outrmse<-c(outrmse,round(rmse(data[,1]+1e-5,data[,2]+1e-5),4))
outslope<-c(outslope,round(as.numeric(linR[[1]][2]),4))
outcc<-c(outcc,round(cor(log10(data[,1]+1e-5),log10(data[,2]+1e-5)),4))


png(paste0(fname,"_prop.png"),width=960,height=960)
par(mar=c(6,6,6,2))
plot(log10(data[,1]),log10(data[,3]),cex.axis=2.5,cex.lab=2.5,cex.main=3,pch=".",xlab=paste0("log10 observed average cut within +-",bias_span,"bp"),ylab=paste0("log10 predicted cut with sum bias in +-",bias_span,"bp / sum bias in +-",tag_span,"bp"),main=paste0("bias window +-",bias_span,"bp\ntotal reads interval +-",tag_span,"bp"))

linR <- lm(log10(data[,3]+1e-5)~log10(data[,1]+1e-5))
abline(linR,lwd=3,col="red")

legend("bottomright",
        c(  paste0("rmse(linear) = ",round(rmse(data[,1]+1e-5,data[,3]+1e-5),4)),
            paste0("slope(log) = ",round(as.numeric(linR[[1]][2]),4)),
            paste0("C.C(log) = ",round(cor(log10(data[,1]+1e-5),log10(data[,3]+1e-5)),4)) ),
        bty="n", cex=3  )
dev.off()
outrmse<-c(outrmse,round(rmse(data[,1]+1e-5,data[,3]+1e-5),4))
outslope<-c(outslope,round(as.numeric(linR[[1]][2]),4))
outcc<-c(outcc,round(cor(log10(data[,1]+1e-5),log10(data[,3]+1e-5)),4))

png(paste0(fname,"_equal.png"),width=960,height=960)
par(mar=c(6,6,6,2))
plot(log10(data[,1]),log10(data[,4]),cex.axis=2.5,cex.lab=2.5,cex.main=3,pch=".",xlab=paste0("log10 observed average cut within +-",bias_span,"bp"),ylab=paste0("log10 predicted cut with equal bias"),main=paste0("bias window +-",bias_span,"bp\ntotal reads interval +-",tag_span,"bp"))

linR <- lm(log10(data[,4]+1e-5)~log10(data[,1]+1e-5))
abline(linR,lwd=3,col="red")

legend("bottomright",
        c(  paste0("rmse(linear) = ",round(rmse(data[,1]+1e-5,data[,4]+1e-5),4)),
            paste0("slope(log) = ",round(as.numeric(linR[[1]][2]),4)),
            paste0("C.C(log) = ",round(cor(log10(data[,1]+1e-5),log10(data[,4]+1e-5)),4)) ),
        bty="n", cex=3  )
dev.off()
outrmse<-c(outrmse,round(rmse(data[,1]+1e-5,data[,4]+1e-5),4))
outslope<-c(outslope,round(as.numeric(linR[[1]][2]),4))
outcc<-c(outcc,round(cor(log10(data[,1]+1e-5),log10(data[,4]+1e-5)),4))

outtable<-rbind(outrmse,outslope,outcc)
colnames(outtable) <- c("pmean","pporp","equal")
write.table(outtable,file=paste0(fname,"_summary.txt"),quote=F,sep="\t",row.names=T,col.names=T)