a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'.pdf')
sudo <- as.numeric(a[3])
data<-read.table(a[1])
each<-(ncol(data)-6)/10
#print(each)
TC <-  apply(data[,(6+1):(6+2*each)],1,mean)
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))
## cleavage
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
## median predict
pl <- data[order(TC),(6+8*each+1):(6+9*each)]
nl <- data[order(TC),(6+9*each+1):(6+10*each)]


pl[(pl<0)] <- 0
nl[(nl<0)] <- 0

pcpl <- log2((pc+sudo) / (pl+sudo))
ncnl <- log2((nc+sudo) / (nl+sudo))
pcpu <- log2((pc+sudo) / (pu+sudo))
ncnu <- log2((nc+sudo) / (nu+sudo))
pcpp <- log2((pc+sudo) / (pp+sudo))
ncnp <- log2((nc+sudo) / (np+sudo))
 ### sitepro
number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)


mpc <- apply(pc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnc <- apply(nc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpb <- apply(pb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnb <- apply(nb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpl <- apply(pl,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnl <- apply(nl,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpu <- apply(pu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnu <- apply(nu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpp <- apply(pp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnp <- apply(np,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)


mpcpl <- apply(pcpl,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnl <- apply(ncnl,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpcpu <- apply(pcpu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnu <- apply(ncnu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpcpp <- apply(pcpp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnp <- apply(ncnp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)


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
    legend("topleft",LE,col=usecolor,lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
}
sitepro2 <- function(md1,md2,n1,n2,usecolor,LE,M){

    plot(md1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M,axes=F,ylim=c(min(md1,md2,n1,n2),max(md2,md1,n1,n2)))
#    print(min(d1,d2))
#    print(max(d1,d2))
    lines(md2,col=usecolor[2])
    lines(n1,col=usecolor[3])
    lines(n2,col=usecolor[4])
    legend("topleft",LE,col=usecolor,lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
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
zmax=max(p80,1)
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



pdf(file=outfile,height=7,width=16)
par(oma = c(1, 1,1, 1),mar=c(1.0, 1.0, 3.0, 1.0))
#layout(matrix(c(1,1,2,2,3,3,4,4,5,7,9,11,13,15,17,19,6,8,10,12,14,16,18,20),nrow=3,byrow=T),heights=c(1,2,0.25))
#layout(matrix(c(1,1,2,2,3,3,4,6,8,10,12,14,5,7,9,11,13,15),nrow=3,byrow=T),heights=c(1,2,0.25))
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,rep(6,5),rep(7,5)),nrow=2,byrow=T),heights=c(1,2))

sitepro(mpc,mnc,c("red","blue"),c("plus","minus"),paste0("Cleavage P : ",round(c_score,4)))
sitepro(mpb,mnb,c("red","blue"),c("plusBG","minusBG"),"Seqbias")
sitepro(mpu,mnu,c("red","blue"),c("plusUniform","minusUniform"),"Uniform")
sitepro(mpp,mnp,c("red","blue"),c("plusProportion","minusProportion"),"Proportion predict")
sitepro(mpl,mnl,c("red","blue"),c("plusMed","minusMed"),"Median predict")

sitepro2(mpcpu - median(c(mpcpu,mncnu)), mncnu- median(c(mpcpu,mncnu)), mpcpl- median(c(mpcpl,mncnl)), mncnl- median(c(mpcpl,mncnl)),c("red","blue","#CDAD00","#228B22"),c("+cut/uniform","-cut/uniform","+cut/median predict","-cut/median predict"),"uniform vs median predict")
sitepro2(mpcpu- median(c(mpcpu,mncnu)), mncnu- median(c(mpcpu,mncnu)), mpcpp- median(c(mpcpp,mncnp)), mncnp- median(c(mpcpp,mncnp)),c("red","blue","#CDAD00","#228B22"),c("+cut/uniform","-cut/uniform","+cut/proportion predict","-cut/proportion predict"),"uniform vs proportion predict")
#sitepro2(mpcpl, mncnl, mpcpp, mncnp,c("red","blue","#CDAD00","#228B22"),c("+cut/uniform","-cut/uniform","+cut/proportion predict","-cut/proportion predict"),"uniform vs proportion predict")


dev.off()
print(c(a[2],c_score))


