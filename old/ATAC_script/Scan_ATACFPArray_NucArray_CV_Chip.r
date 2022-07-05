a<-commandArgs(T)
infile <- a[1]# "~/Desktop/ATAC/Result/FP_nucOC_heatmap/test/CFOS_MA_signal.bed"
outfile <- paste0(a[2],'_cv_chipbw_10k.pdf')

d<-read.table(infile)
if (nrow(d) > 10000){
data <- d[order(d[5],decreasing=T),][1:10000,]
}else{
data<-d
}
print(data[1,1:6])


each<-(ncol(data)-6-400)/6
#print(each)
TC <-  apply(data[,(406+5*each+1):(406+6*each)],1,mean)
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))
atac <- data[order(TC),(6+1):(206)]
mnase <- data[order(TC),(206+1):(406)]
pc <- data[order(TC),(406+1):(406+each)]
nc <- data[order(TC),(406+each+1):(406+each+each)]
pb <- data[order(TC),(406+each+each+1):(406+each+each+each)]
nb <- data[order(TC),(406+3*each+1):(406+4*each)]
cv <- data[order(TC),(406+4*each+1):(406+5*each)]
chip <- data[order(TC),(406+5*each+1):(406+6*each)]

number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)

mpc <- apply(pc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnc <- apply(nc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpb <- apply(pb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnb <- apply(nb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mcv <- apply(cv,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mmnase <- apply(mnase,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
matac <- apply(atac,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mchip <- apply(chip,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
#bottommpc <- apply(pc[1:1000,],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
#bottommnc <- apply(nc[1:1000,],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
#topmpc <- apply(pc[(nrow(pc)-1000):nrow(pc),],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
#topmnc <- apply(nc[(nrow(pc)-1000):nrow(pc),],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
#topmcv <- apply(cv[(nrow(cv)-1000):nrow(cv),],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)
#bottommcv <- apply(cv[1:1000,],2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=10,b=1000)

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
sitepro<-function(md1,md2,usecolor,LE,M){
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
  #  legend("topright",Right,col="black",bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
#    par(new=T)
#    plot(mcv,type="l",col="black",xlab='',ylab="",axes=F)
#    axis(side=4)
    box()
}
sitepro_regular<-function(siteData,usecolor,LE,M){
        #number <- nrow(d1)
        #rm_number <- ceiling(nrow(d1)/200)
        #print(number)
        #print(number-rm_number)
        #md1<-apply(d1,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
        #md2<-apply(d2,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)

    plot(siteData,type="l",col=usecolor,xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=T)
#    print(min(d1,d2))
#    print(max(d1,d2))
#    lines(md2,col=usecolor[2])

    legend("topleft",c(LE),col=c(usecolor),lwd=4,bty="n")
#    legend("topright",Right,col="black",bty="n")
    #xlan <- floor(length(md1)/2)*2
    #axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
#    axis(side=2)
#    par(new=T)
#    plot(mcv,type="l",col="black",xlab='',ylab="",axes=F)
#    axis(side=4)
#    box()
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


bi_heatmap_mnase<-function(data0,usecolor,p,M){
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
#xlan <- floor(ncol(data0)/2)*2
axis(side=1,at = seq(0,2000,500),labels=seq(-1000,1000,500))

box()

#abline(h=29782,lwd=2)
#abline(h=nrow(data0) - 15384 , lwd=2)

image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)

axis(1,seq(round(zmin),round(zmax),round(zmax-zmin)/4),seq(round(zmin),round(zmax),round(zmax-zmin)/4))

box()
}


pdf(file=outfile,height=8.5,width=17)
par(oma = c(1, 1,1, 1),mar=c(1.0, 2, 1.0, 2))
layout(rbind(c(1,1,2,2,3,3,4,4,5,6),c(7,7,9,9,11,13,15,17,19,21),c(8,8,10,10,12,14,16,18,20,22)),heights=c(1,2,0.3))

sitepro_regular(matac,"red","ATAC nuc","")
sitepro_regular(mmnase,"red","MNase","")
sitepro(mpc,mnc,c("red","blue"),c("plus","minus"),paste0("ATAC Cleavage P : ",round(c_score,4)))
sitepro(mpb,mnb,c("red","blue"),c("plusBG","minusBG"),"ATAC Seqbias")
sitepro_regular(mchip,"red","Chip","")
sitepro_regular(mcv,"black","conservation","")


#par(mar=c(0,2,1.0,2))
#sitepro_no(topmpc,topmnc,topmcv,c("red","blue"),c("plus","minus"),"top1k/bottom1k","top1k")
#par(mar=c(1.0,2,0,2))
#sitepro(bottommpc,bottommnc,bottommcv,c("red","blue"),c("plus","minus"),"","bottom1k")

par(mar=c(1.0,2,3,2))
bi_heatmap_mnase(atac,c("white","red"),0.99,"ATAC nuc")
bi_heatmap_mnase(mnase,c("white","red"),0.99,"MNase")
bi_heatmap(pc,c("white","red"),0.99,"ATAC +cleavage",25)
bi_heatmap(nc,c("white","blue"),0.99,"ATAC -cleavage",25)
bi_heatmap(pb,c("white","red"),0.99,"ATAC +Seqbias",25)
bi_heatmap(nb,c("white","blue"),0.99,"ATAC -Seqbias",25)
bi_heatmap(chip,c("white","red"),0.99,"Chipseq",100)
bi_heatmap(cv,c("white","black"),0.99,"Conserv",25)
#par(mfrow=c(1,2))
#bi_heatmap(msuse,c("white","blue"),0.99,"")

dev.off()
#print(c(a[2],c_score))


