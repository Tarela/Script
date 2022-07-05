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
rmse <- mean((as.numeric(data[,1])-as.numeric(data[,3]))**2)**0.5  
sdsummary<-rbind(
make_summary(data[,1]),
make_summary(data[which(data[,3] < 1 &  data[,3] >=0),1]),
make_summary(data[which(data[,3] < 2 &  data[,3] >=1 ),1]),
make_summary(data[which(data[,3] < 3 &  data[,3] >=2 ),1]),
make_summary(data[which(data[,3] < 4 &  data[,3] >=3 ),1]),
make_summary(data[which(data[,3] < 5 &  data[,3] >=4 ),1]),
make_summary(data[which(data[,3] < 6 &  data[,3] >=5 ),1]),
make_summary(data[which(data[,3] < 7 &  data[,3] >=6 ),1]),
make_summary(data[which(data[,3] < 8 &  data[,3] >=7 ),1]),
make_summary(data[which(data[,3] < 9 &  data[,3] >=8 ),1]),
make_summary(data[which(data[,3] < 10 & data[,3] >=9 ),1]),
make_summary(data[which(data[,3] < 11 & data[,3] >=10),1]),
make_summary(data[which(data[,3] < 12 & data[,3] >=11),1]),
make_summary(data[which(data[,3] < 13 & data[,3] >=12),1]),
make_summary(data[which(data[,3] < 14 & data[,3] >=13),1]),
make_summary(data[which(data[,3] < 15 & data[,3] >=14),1]),
make_summary(data[which(data[,3] < 16 & data[,3] >=15),1]),
make_summary(data[which(data[,3] < 17 & data[,3] >=16),1]),
make_summary(data[which(data[,3] < 18 & data[,3] >=17),1]),
make_summary(data[which(data[,3] < 19 & data[,3] >=18),1]),
make_summary(data[which(data[,3] < 20 & data[,3] >=19),1]),
rep(rmse,4))

colnames(sdsummary)<-c("mean","sd","median","mad")
rownames(sdsummary)<-c("all",paste0("pred",seq(0,19,1)),"rmse")
write.table(sdsummary,file=sdtable,quote=F,sep="\t",col.names=T,row.names=T)



pdf(outplot,height=6.5,width=14)
par(mfrow=c(2,2))

i=2
boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),pch=".",ylim=c(0,100),outline=F,xlab="predicted cut",ylab="observed cut",main="original predict")
#box_summary<-boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),plot=F)

#pred_cut_raw <- as.numeric(box_summary$names)[1:30]
#p25_raw <- box_summary$stats[2,1:30]
#p50_raw <- box_summary$stats[3,1:30]
#p75_raw <- box_summary$stats[4,1:30]


i=3
boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),pch=".",ylim=c(0,100),outline=F,xlab="predicted cut",ylab="observed cut",main="forword X rev")


i=4
boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),pch=".",ylim=c(0,100),outline=F,xlab="predicted cut",ylab="observed cut",main="encode log")

i=5
boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),pch=".",ylim=c(0,100),outline=F,xlab="predicted cut",ylab="observed cut",main="encode linear")

#box_summary<-boxplot(data[which(data[,i]<50),1]~floor(data[which(data[,i]<50),i]),plot=F)
#pred_cut_x <- as.numeric(box_summary$names)[1:30]
#p25_x <- box_summary$stats[2,1:30]
#p50_x <- box_summary$stats[3,1:30]
#p75_x <- box_summary$stats[4,1:30]

#summary_table <- rbind(pred_cut_raw,p25_raw,p50_raw,p75_raw,pred_cut_x,p25_x,p50_x,p75_x)
#write.table(summary_table,file=outtable,quote=F,sep="\t",col.names=F,row.names=T)

dev.off()

