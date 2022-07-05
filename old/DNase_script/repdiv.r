a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'_cv_chipbw_10k.pdf')

d<-read.table(infile)
#if (nrow(d) > 10000){
#data <- d[order(d[5],decreasing=T),][1:10000,]
#}else{
data<-d
#}


each<-(ncol(data)-7)/8
TC <-  data[,7]
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))

highoc<- order(TC,decreasing=T)[1:floor(length(TC)/5)]
lowoc<- order(TC,decreasing=F)[1:floor(length(TC)/5)]


hpcC <- as.matrix(data[highoc,(7+1):(7+each)])+0.05
hncC <- as.matrix(data[highoc,(7+each+1):(7+each+each)])+0.05
hpc <- as.matrix(data[highoc,(7+each+each+1):(7+each+each+each)])+0.05
hnc <- as.matrix(data[highoc,(7+3*each+1):(7+4*each)])+0.05
hpc2 <- as.matrix(data[highoc,(7+4*each+1):(7+5*each)])+0.05
hnc2 <- as.matrix(data[highoc,(7+5*each+1):(7+6*each)])+0.05
hpp <- as.matrix(data[highoc,(7+6*each+1):(7+7*each)])+0.05
hnp <- as.matrix(data[highoc,(7+7*each+1):(7+8*each)])+0.05

lpcC <- as.matrix(data[lowoc,(7+1):(7+each)])+0.05
lncC <- as.matrix(data[lowoc,(7+each+1):(7+each+each)])+0.05
lpc <- as.matrix(data[lowoc,(7+each+each+1):(7+each+each+each)])+0.05
lnc <- as.matrix(data[lowoc,(7+3*each+1):(7+4*each)])+0.05
lpc2 <- as.matrix(data[lowoc,(7+4*each+1):(7+5*each)])+0.05
lnc2 <- as.matrix(data[lowoc,(7+5*each+1):(7+6*each)])+0.05
lpp <- as.matrix(data[lowoc,(7+6*each+1):(7+7*each)])+0.05
lnp <- as.matrix(data[lowoc,(7+7*each+1):(7+8*each)])+0.05



hacC <- hpcC+hncC
hac <- hpc+hnc
hac2 <- hpc2+hnc2
hap <- hpp + hnp
lacC <- lpcC + lncC
lac <- lpc+lnc
lac2 <- lpc2+lnc2
lap <- lpp + lnp


#high_cVp <- mean(abs(log2(cbind(hpc,hnc)/cbind(hpp,hnp))))
#high_cVc <- mean(abs(log2(cbind(hpc,hnc)/cbind(hpc2,hnc2))))
#low_cVp <- mean(abs(log2(cbind(lpc,lnc)/cbind(lpp,lnp))))
#low_cVc <- mean(abs(log2(cbind(lpc,lnc)/cbind(lpc2,lnc2))))

high_cVp_all <- mean(abs(log2(hacC/hap)))
high_cVc_all <- mean(abs(log2(hac/hac2)))
low_cVp_all <- mean(abs(log2(lacC/lap)))
low_cVc_all <- mean(abs(log2(lac/lac2)))



#outresult <-c(high_cVp, low_cVp, high_cVc, low_cVc)
outresult <- c(infile,high_cVp_all, low_cVp_all, high_cVc_all, low_cVc_all)
print(outresult)

