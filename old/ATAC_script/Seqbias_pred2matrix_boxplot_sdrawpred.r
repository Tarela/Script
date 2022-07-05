a<-commandArgs(T)
infile <- a[1]
outplot <- paste0(a[2],'.pdf')
#outtable <- paste0(a[2],'_quantile.txt')
sdtable<-paste0(a[2],'_sdmad.txt')
data<-read.table(infile)

make_summary<-function(inputdata){
	outresult <- c(mean(inputdata),sd(inputdata),median(inputdata),mad(inputdata))
	return(outresult)
}
sdsummary<-rbind(
make_summary(data[,1]),
make_summary(data[which(data[,2] < 1 &  data[,2] >=0),1]),
make_summary(data[which(data[,2] < 2 &  data[,2] >=1 ),1]),
make_summary(data[which(data[,2] < 3 &  data[,2] >=2 ),1]),
make_summary(data[which(data[,2] < 4 &  data[,2] >=3 ),1]),
make_summary(data[which(data[,2] < 5 &  data[,2] >=4 ),1]),
make_summary(data[which(data[,2] < 6 &  data[,2] >=5 ),1]),
make_summary(data[which(data[,2] < 7 &  data[,2] >=6 ),1]),
make_summary(data[which(data[,2] < 8 &  data[,2] >=7 ),1]),
make_summary(data[which(data[,2] < 9 &  data[,2] >=8 ),1]),
make_summary(data[which(data[,2] < 10 & data[,2] >=9 ),1]),
make_summary(data[which(data[,2] < 11 & data[,2] >=10),1]),
make_summary(data[which(data[,2] < 12 & data[,2] >=11),1]),
make_summary(data[which(data[,2] < 13 & data[,2] >=12),1]),
make_summary(data[which(data[,2] < 14 & data[,2] >=13),1]),
make_summary(data[which(data[,2] < 15 & data[,2] >=14),1]),
make_summary(data[which(data[,2] < 16 & data[,2] >=15),1]),
make_summary(data[which(data[,2] < 17 & data[,2] >=16),1]),
make_summary(data[which(data[,2] < 18 & data[,2] >=17),1]),
make_summary(data[which(data[,2] < 19 & data[,2] >=18),1]),
make_summary(data[which(data[,2] < 20 & data[,2] >=19),1]))
colnames(sdsummary)<-c("mean","sd","median","mad")
rownames(sdsummary)<-c("all",paste0("pred",seq(0,19,1)))
write.table(sdsummary,file=sdtable,quote=F,sep="\t",col.names=T,row.names=T)



pdf(outplot,height=6.5,width=14)
par(mfrow=c(1,2))

i=2
boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),pch=".",ylim=c(0,100),outline=F,xlab="predicted cut",ylab="observed cut",main="original predict")
#box_summary<-boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),plot=F)

#pred_cut_raw <- as.numeric(box_summary$names)[1:30]
#p25_raw <- box_summary$stats[2,1:30]
#p50_raw <- box_summary$stats[3,1:30]
#p75_raw <- box_summary$stats[4,1:30]


i=3
boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),pch=".",ylim=c(0,100),outline=F,xlab="predicted cut",ylab="observed cut",main="forword X rev")

#box_summary<-boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),plot=F)
#pred_cut_x <- as.numeric(box_summary$names)[1:30]
#p25_x <- box_summary$stats[2,1:30]
#p50_x <- box_summary$stats[3,1:30]
#p75_x <- box_summary$stats[4,1:30]

#summary_table <- rbind(pred_cut_raw,p25_raw,p50_raw,p75_raw,pred_cut_x,p25_x,p50_x,p75_x)
#write.table(summary_table,file=outtable,quote=F,sep="\t",col.names=F,row.names=T)

dev.off()

