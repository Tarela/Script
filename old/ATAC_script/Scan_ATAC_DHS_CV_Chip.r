a<-commandArgs(T)
infile <- a[1]#"~/Desktop/ATAC/Result/FP_nucOC_heatmap/test/NRSF_M00256_K562_signal.bed"
outfile <- a[2]

data<-read.table(infile)


each<-(ncol(data)-5)/6
#print(each)
TC <-  apply(data[,(5+4*each+1):(5+5*each)],1,mean)
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))
atac_dhs <- data[order(TC),(5+1):(5+each)]
atac_mono <- data[order(TC),(5+1*each+1):(5+2*each)]
atac_all <- data[order(TC),(5+2*each+1):(5+3*each)]
dnase <- data[order(TC),(5+3*each+1):(5+4*each)]
chip <- data[order(TC),(5+4*each+1):(5+5*each)]
cv <- data[order(TC),(5+5*each+1):(5+6*each)]


number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)

matac_dhs <- apply(atac_dhs,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
matac_mono <- apply(atac_mono,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
matac_all <- apply(atac_all,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mdnase <- apply(dnase,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mchip <- apply(chip,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mcv<- apply(cv,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)


sitepro_regular<-function(siteData,usecolor,LE,M){

    plot(siteData,type="l",col=usecolor,xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F,lwd=2)
    axis(side=2)
    xlan <- floor(length(siteData)/2)*2
    axis(side=1,at = c(0,100,200,300,400),labels=c(-200,-100,0,100,200))

    box()


    legend("topleft",c(LE),col=c(usecolor),lwd=4,bty="n")
}

 ### heatmap
bi_heatmap<-function(data0,usecolor,p,M,ax){
data<-data0
data <- c(as.matrix(data))
data <- sort(data)
min <- data[1]
max <- data[length(data)]
#print(length(data))
temp<-data[round(c(0.010000,0.5,p)*length(data))]
p20<-temp[1]
p50<-temp[2]
p80<-temp[3]
zmin=p20
zmax=p80
ColorRamp <- colorRampPalette(usecolor, bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
data0[data0<zmin] <- zmin
data0[data0>zmax] <- zmax
ColorRamp_ex <- ColorRamp[round( (min(data0)-zmin)*10000/(zmax-zmin) ) : round( (max(data0)-zmin)*10000/(zmax-zmin) )]
image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab="", ylab="",main=M,useRaster=T)
axis(side=2)
xlan <- floor(ncol(data0)/2)*2
axis(side=1,at = c(0,100,200,300,400),labels=c(-200,-100,0,100,200))

box()

#abline(h=29782,lwd=2)
#abline(h=nrow(data0) - 15384 , lwd=2)

image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)

axis(1,seq(round(zmin),round(zmax),round(zmax-zmin)/4),seq(round(zmin),round(zmax),round(zmax-zmin)/4))

box()
}



pdf(file=outfile,height=8.5,width=17)
par(oma = c(1, 1,1, 1),mar=c(1.0, 2, 1.0, 2))
layout(rbind(c(1,2,3,4,5,6),c(7,9,11,13,15,17),c(8,10,12,14,16,18)),heights=c(1,2,0.3))

sitepro_regular(matac_dhs,"darkred","ATAC le100","")
sitepro_regular(matac_mono,"darkred","ATAC le247","")
sitepro_regular(matac_all,"darkred","ATAC all","")
sitepro_regular(mdnase,"darkred","DNase","")
sitepro_regular(mchip,"darkblue","Chip","")
sitepro_regular(mcv,"black","conservation","")


par(mar=c(1.0,2,3,2))
bi_heatmap(atac_dhs,c("white","darkred"),0.99,"ATAC le100",200)
bi_heatmap(atac_mono,c("white","darkred"),0.99,"ATAC le247",200)
bi_heatmap(atac_all,c("white","darkred"),0.99,"ATAC all",200)
bi_heatmap(dnase,c("white","darkred"),0.99,"DNase",200)
bi_heatmap(chip,c("white","darkblue"),0.95,"Chipseq",200)
bi_heatmap(cv,c("white","black"),0.99,"Conserv",200)
#par(mfrow=c(1,2))
#bi_heatmap(msuse,c("white","blue"),0.99,"")

dev.off()
#print(c(a[2],c_score))


