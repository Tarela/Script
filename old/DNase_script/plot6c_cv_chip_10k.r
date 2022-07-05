a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'_cv_chipbw_10k.pdf')

d<-read.table(a[1])
if (nrow(d) > 10000){
data <- d[order(d[5],decreasing=T),][1:10000,]
}else{
data<-d
}
print(data[1,1:6])


each<-(ncol(data)-8)/8
#print(each)
TC <-  apply(data[,(8+1+7*each):(8+8*each)],1,mean)
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
chip<-data[order(TC),(8+7*each+1):(8+8*each)]
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

bottommpc <- apply(pc[1:1000,],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
bottommnc <- apply(nc[1:1000,],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
topmpc <- apply(pc[(nrow(pc)-1000):nrow(pc),],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
topmnc <- apply(nc[(nrow(pc)-1000):nrow(pc),],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
topmcv <- apply(cv[(nrow(cv)-1000):nrow(cv),],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
bottommcv <- apply(cv[1:1000,],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)

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
sitepro<-function(md1,md2,mcv,usecolor,LE,M,Right){
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
    legend("topright",Right,col="black",bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    par(new=T)
    plot(mcv,type="l",col="black",xlab='',ylab="",axes=F)
    axis(side=4)
    box()
}
sitepro_no<-function(md1,md2,mcv,usecolor,LE,M,Right){
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
    legend("topright",Right,col="black",bty="n")
    #xlan <- floor(length(md1)/2)*2
    #axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    par(new=T)
    plot(mcv,type="l",col="black",xlab='',ylab="",axes=F)
    axis(side=4)
    box()
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
axis(side=1,at = seq(0,xlan,25),labels=seq(-ax,ax,ax))

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



pdf(file=outfile,height=9,width=17)
par(oma = c(1, 1,1, 1),mar=c(1.0, 2, 1.0, 2))
layout(matrix(c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,5,5,6,8,10,12,14,16,18,20,7,9,11,13,15,17,19,21),nrow=4,byrow=T),heights=c(0.5,0.5,2,0.25))
sitepro(mpc,mnc,mcv,c("red","blue"),c("plus","minus"),paste0("Cleavage P : ",round(c_score,4)),"")
sitepro(mpi,mni,mcv,c("red","blue"),c("plus","minus"),"Naked DNA","")
sitepro(mpb,mnb,mcv,c("red","blue"),c("plusBG","minusBG"),"Seqbias","")
par(mar=c(0,2,1.0,2))
sitepro_no(topmpc,topmnc,topmcv,c("red","blue"),c("plus","minus"),"top1k/bottom1k","top1k")
par(mar=c(1.0,2,0,2))
sitepro(bottommpc,bottommnc,bottommcv,c("red","blue"),c("plus","minus"),"","bottom1k")

par(mar=c(1.0,2,3,2))
bi_heatmap(pc,c("white","red"),0.99,"+cleavage",25)
bi_heatmap(nc,c("white","blue"),0.99,"-cleavage",25)
bi_heatmapInput(pI,c("white","red"),0.99,"+Naked")
bi_heatmapInput(nI,c("white","blue"),0.99,"-Naked")
bi_heatmap(pb,c("white","red"),0.99,"+Seqbias",25)
bi_heatmap(nb,c("white","blue"),0.99,"-Seqbias",25)
bi_heatmap(chip,c("white","red"),0.99,"Chipseq",100)
bi_heatmap(cv,c("white","black"),0.99,"Conserv",25)
#par(mfrow=c(1,2))
#bi_heatmap(msuse,c("white","blue"),0.99,"")

dev.off()
print(c(a[2],c_score))

