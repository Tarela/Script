a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'_cv.pdf')
data<-read.table(a[1])

each<-(ncol(data)-8)/7
#print(each)
TC <-  apply(data[,(8+1):(8+2*each)],1,mean)
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))
pc <- data[order(TC),(8+1):(8+each)]
nc <- data[order(TC),(8+each+1):(8+each+each)]
pI <- data[order(TC),(8+each+each+1):(8+each+each+each)]
nI <- data[order(TC),(8+3*each+1):(8+4*each)]
pb <- data[order(TC),(8+4*each+1):(8+5*each)]
nb <- data[order(TC),(8+5*each+1):(8+6*each)]
cv <- data[order(TC),(8+6*each+1):(8+7*each)]

tagcount <- data[order(TC),7]
FOS<- data[order(TC),8]
ms <- data[order(TC),5]
msuse<-cbind(ms,ms,ms)
#data[order(TC),1:6][1:100,]
 ### sitepro
number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)

mpc <- apply(pc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnc <- apply(nc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpi <- apply(pI,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mni <- apply(nI,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpb <- apply(pb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnb <- apply(nb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mcv <- apply(cv,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
#par(mfrow=c(2,1))
#plot(mpc,type="l",col="red")
#lines(mnc,col="blue")
#plot(mpb,type="l",col="red")
#lines(mnb,col="blue")

cmp_rg <- (floor(ncol(pc)/2)-24 ):(floor(ncol(pc)/2)+25 ) 
fore<-c(mpc[cmp_rg],mnc[cmp_rg])
back<-c(mpb[cmp_rg],mnb[cmp_rg]) 
###make porportion
foreP <- fore/sum(fore)
backP <- back/sum(back)
#c_score <- cor(foreP,backP)
#print(length(mpc))
#print(length(cmp_rg))
c_score<-cor(c(mpc,mnc),c(mpb,mnb))
#print(c_score)

#m_score <- sum(foreP*log(foreP/backP))
#a<-pc[,100]
## sitepor
sitepro<-function(md1,md2,mcv,usecolor,LE,M){
	#number <- nrow(d1)
	#rm_number <- ceiling(nrow(d1)/200)
	#print(number)
	#print(number-rm_number)
	#md1<-apply(d1,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
	#md2<-apply(d2,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
	
    plot(md1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F,ylim=c(min(md1,md2),max(md2,md1)))
#    print(min(d1,d2))
#    print(max(d1,d2))
    lines(md2,col=usecolor[2])

    legend("topleft",c(LE,"conserv"),col=c(usecolor,"black"),lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    par(new=T)
    plot(mcv,type="l",col="black",xlab='',ylab="",axes=F)
    axis(side=4)
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

axis(1,seq(round(zmin),round(zmax),round(zmax-zmin)/4),seq(round(zmin),round(zmax),round(zmax-zmin)/4))

box()
}
bi_heatmapInput<-function(data0,usecolor,p,M){
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
zmin=0#p20
zmax=1#p80
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
axis(1,c(0,1),c(0,1))
box()
}



pdf(file=outfile,height=9,width=16)
par(oma = c(1, 1,1, 1),mar=c(1.0, 1.0, 1.0, 1.0))
layout(matrix(c(1,1,2,2,3,3,4,6,8,10,12,14,5,7,9,11,13,15),nrow=3,byrow=T),heights=c(1,2,0.25))
sitepro(mpc,mnc,mcv,c("red","blue"),c("plus","minus"),paste0("Cleavage P : ",round(c_score,4)))
sitepro(mpi,mni,mcv,c("red","blue"),c("plus","minus"),"Naked DNA")
sitepro(mpb,mnb,mcv,c("red","blue"),c("plusBG","minusBG"),paste0("average FP score : ",round(mean(FOS),4)))

bi_heatmap(pc,c("white","red"),0.99,"")
bi_heatmap(nc,c("white","blue"),0.99,"")
bi_heatmapInput(pI,c("white","red"),0.99,"")
bi_heatmapInput(nI,c("white","blue"),0.99,"")
bi_heatmap(pb,c("white","red"),0.99,"")
bi_heatmap(nb,c("white","blue"),0.99,"")
#par(mfrow=c(1,2))
#bi_heatmap(msuse,c("white","blue"),0.99,"")

dev.off()
print(c(a[2],c_score))

