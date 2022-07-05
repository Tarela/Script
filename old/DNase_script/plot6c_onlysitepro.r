a<-commandArgs(T)
infile <- a[1]
outfile <- paste0(a[2],'_onlysitepro.pdf')
data<-read.table(a[1])

each<-(ncol(data)-8)/6
print(each)
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
pb[pb == 0] <- 0.001
nb[nb == 0] <- 0.001
psub <- pc/pb 
nsub <- nc/nb
print(c(min(pc),min(pI),min(pb)))
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
mpsub <- apply(psub,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnsub <- apply(nsub,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)

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
print(length(mpc))
print(length(cmp_rg))
c_score<-cor(c(mpc,mnc),c(mpb,mnb))
print(c_score)

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

pdf(file=outfile,height=3,width=16)
par(oma = c(1, 1,1, 1),mar=c(1.0, 1.0, 1.0, 1.0),mfrow=c(1,4))
#layout(matrix(c(1,1,2,2,3,3,4,4,5,7,9,11,13,15,17,19,6,8,10,12,14,16,18,20),nrow=3,byrow=T),heights=c(1,2,0.25))

sitepro(mpc,mnc,c("red","blue"),c("plus","minus"),paste0("Cleavage P : ",round(c_score,4)))
sitepro(mpi,mni,c("red","blue"),c("plus","minus"),"Naked DNA")
sitepro(mpb,mnb,c("red","blue"),c("plusBG","minusBG"),"Seqbias")

sitepro(mpsub,mnsub,c("red","blue"),c("plusCut/plusBG","minusCut/minusBG"),"Cut/BG")


dev.off()
print(c(a[2],c_score))

