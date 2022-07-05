a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'.pdf')
sudo <- 1
d<-read.table(infile)
if (nrow(d) > 10000){
data <- d[order(d[5],decreasing=T),][1:10000,]
}else{
data<-d
}
#print(data[1,1:6])


each<-(ncol(data)-6)/10
#print(each)
TC <-  apply(data[,(6+1+9*each):(6+10*each)],1,mean)
pc <- data[order(TC),(6+1):(6+each)]
nc <- data[order(TC),(6+each+1):(6+each+each)]
## bias
pb <- data[order(TC),(6+each+each+1):(6+each+each+each)]
nb <- data[order(TC),(6+3*each+1):(6+4*each)]
## uniform
pu <- data[order(TC),(6+4*each+1):(6+5*each)]
nu <- data[order(TC),(6+5*each+1):(6+6*each)]
## proportion predict
pp <- data[order(TC),(6+6*each+1):(6+7*each)]
np <- data[order(TC),(6+7*each+1):(6+8*each)]
## cv and chip
cv <- data[order(TC),(6+8*each+1):(6+9*each)]
chip <- data[order(TC),(6+9*each+1):(6+10*each)]

pcpu <- log((pc+sudo) / (pu+sudo))
ncnu <- log((nc+sudo) / (nu+sudo))
pcpp <- log((pc+sudo) / (pp+sudo))
ncnp <- log((nc+sudo) / (np+sudo))

 ### sitepro
number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)

mpcpu <- apply(pcpu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnu <- apply(ncnu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpcpp <- apply(pcpp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnp <- apply(ncnp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mcv <- apply(cv,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)


## sitepor
sitepro0<-function(md1,usecolor,LE,M){
	#number <- nrow(d1)
	#rm_number <- ceiling(nrow(d1)/200)
	#print(number)
	#print(number-rm_number)
	#md1<-apply(d1,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
	#md2<-apply(d2,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
	
    plot(md1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F)
#    print(min(d1,d2))
#    print(max(d1,d2))
    legend("topleft",LE,col=usecolor,lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
}
sitepro2 <- function(md1,md2,n1,n2,usecolor,LE,M){

    plot(md1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M[1],axes=F,ylim=c(min(md1,md2,n1,n2),max(md2,md1,n1,n2)))
#    print(min(d1,d2))
#    print(max(d1,d2))
    lines(md2,col=usecolor[2])
    legend("topleft",LE[1:2],col=usecolor[1:2],lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
    plot(n1,type="l",col=usecolor[3],xlab="Relative distance (bp)",ylab="Average profile",main=M[2],axes=F,ylim=c(min(md1,md2,n1,n2),max(md2,md1,n1,n2)))
    lines(n2,col=usecolor[4])
    legend("topleft",LE[3:4],col=usecolor[3:4],lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
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



pdf(file=outfile,height=9,width=17)
par(oma = c(1, 1,1, 1),mar=c(1.0, 2, 1.0, 2))
layout(matrix(c(1,1,2,2,3,3,4,6,8,10,12,14,5,7,9,11,13,15),nrow=3,byrow=T),heights=c(1,2,0.25))
sitepro2(mpcpu - median(c(mpcpu,mncnu)), mncnu- median(c(mpcpu,mncnu)), mpcpp- median(c(mpcpp,mncnp)), mncnp- median(c(mpcpp,mncnp)),c("red","blue","red","blue"),c("+cut/uniform","-cut/uniform","+cut/proportion bias","-cut/proportion bias"),c("uniform","proportion"))
sitepro0(mcv,c("black"),c("conservation"),"conservation")

par(mar=c(1.0,2,3,2))
bi_heatmap(pcpu,c("white","red"),0.99,"+cut/uniform",25)
bi_heatmap(ncnu,c("white","blue"),0.99,"-cut/uniform",25)
bi_heatmap(pcpp,c("white","red"),0.99,"+cut/proportion",25)
bi_heatmap(ncnp,c("white","blue"),0.99,"-cut/proportion",25)
bi_heatmap(chip,c("white","red"),0.99,"Chipseq",100)
bi_heatmap(cv,c("white","black"),0.99,"Conserv",25)
#par(mfrow=c(1,2))
#bi_heatmap(msuse,c("white","blue"),0.99,"")

dev.off()
