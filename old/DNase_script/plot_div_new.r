a<-commandArgs(T)
infile <- a[1]	
outfile <- paste0(a[2],'_div.pdf')

d<-read.table(infile)
#if (nrow(d) > 10000){
#data <- d[order(d[5],decreasing=T),][1:10000,]
#}else{
data<-d
#}


each<-(ncol(data)-7)/6
TC <-  data[,7]
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))

oc<- order(TC)


pc <- as.matrix(data[oc,(7+1):(7+each)])+0.05
nc <- as.matrix(data[oc,(7+each+1):(7+each+each)])+0.05
allc <- pc+nc
pc2 <- as.matrix(data[oc,(7+each+each+1):(7+each+each+each)])+0.05
nc2 <- as.matrix(data[oc,(7+3*each+1):(7+4*each)])+0.05
allc2 <- pc2+nc2
pp <- as.matrix(data[oc,(7+4*each+1):(7+5*each)])+0.05
np <- as.matrix(data[oc,(7+5*each+1):(7+6*each)])+0.05
allp <- pp+np

allcvp <- abs(log2(allc/allp))
allcvc <- abs(log2(allc/allc2))
ave_allcvp <- cbind(apply(allcvp,1,mean),apply(allcvp,1,mean))
ave_allcvc <- cbind(apply(allcvc,1,mean),apply(allcvc,1,mean))

chip <- cbind(TC[order(TC)],TC[order(TC)])

number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)

mpc <- apply(pc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnc <- apply(nc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpc2 <- apply(pc2,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnc2 <- apply(nc2,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpp <- apply(pp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnp <- apply(np,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mallc <- apply(allc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mallc2 <- apply(allc2,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mallp <- apply(allp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mallcvp <- apply(allcvp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mallcvc <- apply(allcvc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)

sitepro<-function(md1,md2,usecolor,LE,M){
    plot(md1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F,ylim=c(min(md1,md2),max(md2,md1)))
    lines(md2,col=usecolor[2])
    legend("topleft",LE,col=usecolor,lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
#    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
}

sitepro1<-function(md1,usecolor,LE,M){
    plot(md1,type="l",col=usecolor,xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F)
    #lines(md2,col=usecolor[2])
    legend("topleft",LE,col=usecolor,lwd=4,bty="n")
    #xlan <- floor(length(md1)/2)*2
#    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
}

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
#xlan <- floor(ncol(data0)/2)*2
#axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))

box()

#abline(h=29782,lwd=2)
#abline(h=nrow(data0) - 15384 , lwd=2)

image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)

axis(1,seq(round(zmin),round(zmax),round(zmax-zmin)/4),seq(round(zmin),round(zmax),round(zmax-zmin)/4))

box()
}


pdf(file=outfile,height=6,width=9)
par(oma = c(1, 1,1, 1),mar=c(1.0, 1.0, 1.0, 1.0))
layout(matrix(c(1,2,16,3,17,18,4,6,8,10,12,14,5,7,9,11,13,15),nrow=3,byrow=T),heights=c(1,2,0.25),width=c(1,1,0.2,1,0.2,0.2))
sitepro1(mallc,"red","","observe cleavage,mean")
sitepro1(mallcvp,"blue","","abs(log2(obs/exp)),mean")
sitepro1(mallcvc,"blue","","abs(log2(rep1/rep2)),mean")

bi_heatmap(allc,c("white","red"),0.99,"")
bi_heatmap(allcvp,c("white","blue"),0.99,"")
bi_heatmap(ave_allcvp,c("white","blue"),0.99,"")
bi_heatmap(allcvc,c("white","blue"),0.99,"")
bi_heatmap(ave_allcvc,c("white","blue"),0.99,"")

bi_heatmap(chip,c("white","brown"),0.99,"Chip")
plot(1,1,type="n",axes=F,xlab="",ylab="")
plot(1,1,type="n",axes=F,xlab="",ylab="")
plot(1,1,type="n",axes=F,xlab="",ylab="")
dev.off()
