a<-commandArgs(T)
infile <- a[1]#"~/Desktop/ATAC/Result/FP_nucOC_heatmap/test/NRSF_M00256_K562_signal.bed"
outfile <- paste0(a[2],'.pdf')

d<-read.table(infile)
if (nrow(d) > 10000){
data <- d[order(d[5],decreasing=T),][1:10000,]
}else{
data<-d
}
#print(data[1,1:6])


each<-(ncol(data)-6)/11
#print(each)
TC <-  apply(data[,(6+2*each+1):(6+4*each)],1,mean)
#length(which(TC > 0.5))
#hist(TC*50,nclass=200,xlim=c(0,50))
#plot(density(TC))
pac <- data[order(TC),(6+1):(6+each)][,1:92]
nac <- data[order(TC),(6+each+1):(6+each*2)][,9:100]
#pas <- data[order(TC),(6+2*each+1):(6+each*3)]
#nas <- data[order(TC),(6+3*each+1):(6+4*each)]
#pab <- data[order(TC),(6+4*each+1):(6+5*each)] 
#nab  <- data[order(TC),(6+5*each+1):(6+6*each)]
pdc <- data[order(TC),(6+6*each+1):(6+7*each)]
ndc <- data[order(TC),(6+7*each+1):(6+8*each)]
pdb <- data[order(TC),(6+8*each+1):(6+9*each)]
ndb <- data[order(TC),(6+9*each+1):(6+10*each)]
cv <- data[order(TC),(6+10*each+1):(6+11*each)]

number <- nrow(data)
rm_number <- ceiling(nrow(data)/100)

mpac <- apply(pac,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mnac <- apply(nac,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
#mpas <- apply(pas,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
#mnas <- apply(nas,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
#mpab <- apply(pab,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
#mnab <- apply(nab,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpdc <- apply(pdc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mndc <- apply(ndc,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mpdb <- apply(pdb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mndb <- apply(ndb,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)
mcv  <- apply(cv,2,function(x,a,b) mean(x[order(x)[a:(b-a)]]),a=rm_number,b=number)



#par(mfrow=c(2,1))
#plot(mpc,type="l",col="red")
#lines(mnc,col="blue")
#plot(mpb,type="l",col="red")
#lines(mnb,col="blue")

cmp_rg <- (floor(ncol(pdc)/2)-49 ):(floor(ncol(pdc)/2)+50 ) 
#foreA<-c(mpac[cmp_rg],mnac[cmp_rg])
#foreA_short<-c(mpas[cmp_rg],mnas[cmp_rg])
#backA<-c(mpab[cmp_rg],mnab[cmp_rg]) 
foreD <- c(mpdc[cmp_rg],mndc[cmp_rg])
backD<-c(mpdb[cmp_rg],mndb[cmp_rg]) 

###make porportion
#foreAP <- foreA/sum(foreA)
#foreAP_short <- foreA_short/sum(foreA_short)
#backAP <- backA/sum(backA)
foreDP <- foreD/sum(foreD)
backDP <- backD/sum(backD)

#c_score <- cor(foreP,backP)
#print(length(mpc))
#print(length(cmp_rg))
#c_AP<-cor(foreAP,backAP)
#c_AP_short <- cor(foreAP_short,backAP)
c_DP <- cor(foreDP,backDP)
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

pdf(file=outfile,height=8.5,width=17)
par(oma = c(1, 1,1, 1),mar=c(1.0, 2, 1.0, 2))
layout(rbind(c(1,1,2,2,3,3,4),seq(5,17,2),seq(6,18,2)),heights=c(1,2,0.3))


sitepro(mpac, mnac,c("red","blue"),c("+ATAC","-ATAC"),"ATAC")
#sitepro(mpas, mnas,c("red","blue"),c("+ATAC<100","-ATAC<100"),paste0("ATAC<100 : ",round(c_AP_short,4)))
#sitepro(mpab, mnab,c("red","blue"),c("+ATACbg","-ATACbg"),"ATAC bias")
sitepro(mpdc, mndc,c("red","blue"),c("+DNase","-DNase"),"DNase")
sitepro(mpdb, mndb,c("red","blue"),c("+DNasebg","-DNasebg"),"DNase bias")
sitepro_regular(mcv,"black","conservation","")

par(mar=c(1.0,2,3,2))
bi_heatmap(pac,c("white","red"),0.99,"+ATAC")
bi_heatmap(nac,c("white","blue"),0.99,"-ATAC")
#bi_heatmap(pas,c("white","red"),0.99,"+ATAC<100")
#bi_heatmap(nas,c("white","blue"),0.99,"-ATAC<100")
#bi_heatmap(pab,c("white","red"),0.99,"+ATACbg")
#bi_heatmap(nab,c("white","blue"),0.99,"-ATACbg")
bi_heatmap(pdc,c("white","red"),0.99,"+DNase")
bi_heatmap(ndc,c("white","blue"),0.99,"-DNase")
bi_heatmap(pdb,c("white","red"),0.99,"+DNasebg")
bi_heatmap(ndb,c("white","blue"),0.99,"-DNasebg")
bi_heatmap(cv,c("white","black"),0.99,"Conserv")
#par(mfrow=c(1,2))
#bi_heatmap(msuse,c("white","blue"),0.99,"")

dev.off()
#print(c(a[2],nrow(d), c_AP, c_AP_short, c_DP))


