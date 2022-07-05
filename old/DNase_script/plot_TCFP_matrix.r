a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'.pdf')
data<-read.table(a[1])


each<-(ncol(data)-8)/4
TC <- apply(data[,(8+1):(8+2*each)],1,mean)
pc <- data[order(TC),(8+1):(8+each)]
nc <- data[order(TC),(8+each+1):(8+each+each)]
pb <- data[order(TC),(8+each+each+1):(8+each+each+each)]
nb <- data[order(TC),(8+3*each+1):(8+4*each)]

tagcount<-data[order(TC),7]
FOS<-data[order(TC),8]

 ### sitepro
mpc <- apply(pc,2,mean)
mnc <- apply(nc,2,mean)
mpb <- apply(pb,2,mean)
mnb <- apply(nb,2,mean)

#par(mfrow=c(2,1))
#plot(mpc,type="l",col="red")
#lines(mnc,col="blue")
#plot(mpb,type="l",col="red")
#lines(mnb,col="blue")

cmp_rg <- (floor(ncol(pc)/2)-24 ):(floor(ncol(pc)/2)+25 ) 
fore<-c(mpc[cmp_rg],mnc[cmp_rg])
back<-c(mpb[cmp_rg],mnb[cmp_rg]) ###make porportion
foreP <- fore/sum(fore)
backP <- back/sum(back)
### Correlation coefficient
c_score <- cor(foreP,backP)
#m_score <- sum(foreP*log(foreP/backP))

#info<-c(a[2],)

## sitepor
sitepro<-function(d1,d2,usecolor,LE,M){
    plot(d1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F,ylim=c(min(d1,d2),max(d2,d1)))
#    print(min(d1,d2))
#    print(max(d1,d2))
    lines(d2,col=usecolor[2])
    legend("topleft",LE,col=usecolor,lwd=4,bty="n")
    xlan <- floor(length(mpc)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
}

 ### heatmap
bi_heatmap<-function(data0,usecolor,p,M){
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
axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
box()
#abline(h=29782,lwd=2)
#abline(h=nrow(data0) - 15384 , lwd=2)
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)
#print(zmax)
#print(zmin)
#print(round(zmax-zmin)/5)
#print(seq(round(zmin),round(zmax),round(zmax-zmin)/5))
axis(1,seq(round(zmin),round(zmax),round(zmax-zmin,4)/4),seq(round(zmin),round(zmax),round(zmax-zmin,4)/4))
box()
}
pdf(file=outfile)
par(oma = c(1, 1,1, 1),mar=c(1.0, 1.0, 1.0, 1.0))
layout(matrix(c(1,1,2,2,3,5,7,9,4,6,8,10),nrow=3,byrow=T),heights=c(1,2,0.25))
sitepro(mpc,mnc,c("red","blue"),c("plus","minus"),paste0("Pearsonr : ",round(c_score,4)))
sitepro(mpb,mnb,c("red","blue"),c("plusBG","minusBG"),paste0("meanFOS : ",round(mean(FOS),4)))

bi_heatmap(pc,c("white","red"),0.95,"")
bi_heatmap(nc,c("white","blue"),0.95,"")
bi_heatmap(pb,c("white","red"),0.95,"")
bi_heatmap(nb,c("white","blue"),0.95,"")

dev.off()

