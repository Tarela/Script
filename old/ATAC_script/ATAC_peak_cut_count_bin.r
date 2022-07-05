a<-commandArgs(T)


outname<-a[2]
infile<-a[1]

data <- read.table(infile)
data<-data[order(data[,6],decreasing=T),]
data<-data[which(data[,1]!='chrM'),]

l2<-length(data[which(cumsum(data[,6])<= sum(data[,6])*0.2),6])
l4<-length(data[which(cumsum(data[,6])<= sum(data[,6])*0.4 &cumsum(data[,6])> sum(data[,6])*0.2),6])
l6<-length(data[which(cumsum(data[,6])<= sum(data[,6])*0.6 &cumsum(data[,6])> sum(data[,6])*0.4),6])
l8<-length(data[which(cumsum(data[,6])<= sum(data[,6])*0.8 &cumsum(data[,6])> sum(data[,6])*0.6),6])
l10<-length(data[which(cumsum(data[,6])<= sum(data[,6])*1 &cumsum(data[,6])> sum(data[,6])*0.8),6])


q2<-which(cumsum(data[,6])<= sum(data[,6])*0.2)
q4<-which(cumsum(data[,6])<= sum(data[,6])*0.4 &cumsum(data[,6])> sum(data[,6])*0.2)
q6<-which(cumsum(data[,6])<= sum(data[,6])*0.6 &cumsum(data[,6])> sum(data[,6])*0.4)
q8<-which(cumsum(data[,6])<= sum(data[,6])*0.8 &cumsum(data[,6])> sum(data[,6])*0.6)
q10<-which(cumsum(data[,6])<= sum(data[,6])*1.0 &cumsum(data[,6])> sum(data[,6])*0.8)

rain<-rainbow(5)
pdf(file=paste0(outname,"_nochrM_cutcountCDF.pdf"))
plot(ecdf(log10(data[,6])),verticals=T,do.points=F,main="cumulate log10 cut count in all ATAC peaks",ylab="cumlate",xlab="log10 cut count")

abline(v=log10(data[1,6]),col=rain[1])
abline(v=log10(data[q2[length(q2)],6]),col=rain[2])
abline(v=log10(data[q4[length(q4)],6]),col=rain[3])
abline(v=log10(data[q6[length(q6)],6]),col=rain[4])
abline(v=log10(data[q8[length(q8)],6]),col=rain[5])
abline(v=log10(data[q10[length(q10)],6]))

legend("topleft",paste(c("peak number 0~20% = ","peak number 20~40% = ","peak number 40~60% = ","peak number 60~80% = ","peak number 80~100% = "),c(l2,l4,l6,l8,l10)),col=rainbow(5),lwd=3,bty="n")

dev.off()

write.table(data[q2,1:5],file=paste0(outname,"_nochrM_q2.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(data[q4,1:5],file=paste0(outname,"_nochrM_q4.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(data[q6,1:5],file=paste0(outname,"_nochrM_q6.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(data[q8,1:5],file=paste0(outname,"_nochrM_q8.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(data[q10,1:5],file=paste0(outname,"_nochrM_q10.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(data[,1:5],file=paste0(outname,"_nochrM.bed"),sep="\t",quote=F,row.names=F,col.names=F)


