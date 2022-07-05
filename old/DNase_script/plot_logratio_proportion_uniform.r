a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'.pdf')
sudo <- 1#as.numeric(a[3])
data<-read.table(a[1])
each<-(ncol(data)-6)/8
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

## log ratio
pcpu <- log((pc+sudo) / (pu+sudo))
ncnu <- log((nc+sudo) / (nu+sudo))
pcpp <- log((pc+sudo) / (pp+sudo))
ncnp <- log((nc+sudo) / (np+sudo))
 ### sitepro
number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)


mpc <- apply(pc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnc <- apply(nc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpb <- apply(pb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnb <- apply(nb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpu <- apply(pu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnu <- apply(nu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpp <- apply(pp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnp <- apply(np,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)

mpcpu <- apply(pcpu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnu <- apply(ncnu,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpcpp <- apply(pcpp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mncnp <- apply(ncnp,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)


cmp_rg <- (floor(ncol(pc)/2)-24 ):(floor(ncol(pc)/2)+25 ) 
fore<-c(mpc[cmp_rg],mnc[cmp_rg])
back<-c(mpb[cmp_rg],mnb[cmp_rg]) 
###make porportion
foreP <- fore/sum(fore)
backP <- back/sum(back)
c_score<-cor(c(mpc,mnc),c(mpb,mnb))

### varience
logratio <- c(mpcpp[cmp_rg], mncnp[cmp_rg])
LRscore <- sum( (logratio - mean(logratio))**2 )
#print(length(logratio))
#print(LRscore/(length(logratio)-1))
#print(var(logratio))
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
sitepro3 <- function(md1,md2,n1,n2,usecolor,LE,M){

    plot(md1,type="l",col=usecolor[1],xlab="Relative distance (bp)",ylab="Average profile",main=M[1],axes=F,ylim=c(-1,1.5))
#    print(min(d1,d2))
#    print(max(d1,d2))
    lines(md2,col=usecolor[2])
    legend("topleft",LE[1:2],col=usecolor[1:2],lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
    plot(n1,type="l",col=usecolor[3],xlab="Relative distance (bp)",ylab="Average profile",main=M[2],axes=F,ylim=c(-1,1.5))
    lines(n2,col=usecolor[4])
    legend("topleft",LE[3:4],col=usecolor[3:4],lwd=4,bty="n")
    xlan <- floor(length(md1)/2)*2
    axis(side=1,at = seq(0,xlan,25),labels=seq(-xlan/2,xlan/2,25))
    axis(side=2)
    box()
}

pdf(file=outfile,height=7,width=16)
par(oma = c(1, 1,1, 1),mar=c(1.0, 1.0, 3.0, 1.0))
#layout(matrix(c(1,1,2,2,3,3,4,4,5,7,9,11,13,15,17,19,6,8,10,12,14,16,18,20),nrow=3,byrow=T),heights=c(1,2,0.25))
#layout(matrix(c(1,1,2,2,3,3,4,6,8,10,12,14,5,7,9,11,13,15),nrow=3,byrow=T),heights=c(1,2,0.25))
layout(matrix(c(1,2,3,4,rep(5,2),rep(6,2)),nrow=2,byrow=T),heights=c(1,2))

sitepro(mpc,mnc,c("red","blue"),c("plus","minus"),paste0("Cleavage P : ",round(c_score,4)))
sitepro(mpb,mnb,c("red","blue"),c("plusBG","minusBG"),"Seqbias")
sitepro(mpu,mnu,c("red","blue"),c("plusUniform","minusUniform"),"Uniform")
sitepro(mpp,mnp,c("red","blue"),c("plusProportion","minusProportion"),"Proportion predict")

sitepro2(mpcpu - median(c(mpcpu,mncnu)), mncnu- median(c(mpcpu,mncnu)), mpcpp- median(c(mpcpp,mncnp)), mncnp- median(c(mpcpp,mncnp)),c("red","blue","red","blue"),c("+cut/uniform","-cut/uniform","+cut/proportion bias","-cut/proportion bias"),c("uniform",paste0("proportion bias , logratioScore = ",round(LRscore,4))))
#sitepro3(mpcpu - median(c(mpcpu,mncnu)), mncnu- median(c(mpcpu,mncnu)), mpcpp- median(c(mpcpp,mncnp)), mncnp- median(c(mpcpp,mncnp)),c("red","blue","red","blue"),c("+cut/uniform","-cut/uniform","+cut/proportion bias","-cut/proportion bias"),c("uniform",paste0("proportion bias , logratioScore = ",round(LRscore,4))))


dev.off()
print(c(a[2],LRscore))


