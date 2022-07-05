a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- as.numeric(a[2])#1
datatype <- a[3]#"ATAC"
iter = 1000

ATAC_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_",datatype,"_readscount_seqtype.bed"),header=T,row.names=4)
seqtype_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T,row.names=4)

if(datatype == "ATAC"){
    bias <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/prediction_separate_reads_cutseq/MLE_model/sampling/biasMat/YeastNaked_ATACPE_enc",flankN,"flank.txt"),row.names=1,header=T)[,2]
#    simplex <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Data/bias_matrix/simplex_parameters/ATAC_flank",flankN,"_simplexParam.txt"))[,1]
}else{
    bias <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/prediction_separate_reads_cutseq/MLE_model/sampling/biasMat/IMR90naked_DNaseSE_enc",flankN,"flank.txt"),row.names=1,header=T)[,2]
#    simplex <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Data/bias_matrix/simplex_parameters/DNase_flank",flankN,"_simplexParam.txt"))[,1]    
}

disTSS <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_ADmergePeak_summit_disTSS.bed"),row.names=4)
disTTS <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_ADmergePeak_summit_disTTS.bed"),row.names=4)

disTrans <- function(indis){
  outdis <- abs(indis)
  outdis[which(outdis==0)]<-1
  usedis <- log10(outdis)
  #usedis[which(indis < 0)] <- -1 * usedis[which(indis < 0)] 
  return(usedis)
}
useTSS <- disTrans(disTSS[,10])
useTTS <- disTrans(disTTS[,10])
names(useTSS) <- rownames(disTSS)
names(useTTS) <- rownames(disTTS)

TSSpeak <- names(useTSS)[useTSS <= log10(5000) & useTTS >= log10(5000)]
TTSpeak <- names(useTSS)[useTTS <= log10(5000) & useTSS >= log10(5000)]
otherpeak <- names(useTSS)[useTTS >= log10(5000) & useTSS >= log10(5000)]


X <- ATAC_raw[,4:ncol(ATAC_raw)]#[1:1000,]
L <- seqtype_raw[,4:ncol(seqtype_raw)]#[1:1000,]

X = matrix(as.numeric(as.matrix(X)), ncol = ncol(X))
L = matrix(as.numeric(as.matrix(L)), ncol = ncol(L))
rownames(X) <- rownames(ATAC_raw)
rownames(L) <- rownames(ATAC_raw)
pos = which(rowSums(X) > 0 & rowSums(X) < 1000)

X = X[pos, ]
L = L[pos, ]
#H <- matrix(as.numeric(as.matrix(Hmat)), ncol = ncol(Hmat))[,2:ncol(Hmat)]

vector_generator = function(X, L, iter){
  
  Amat_sumi <- apply(X,2,sum)
  Amat_sumj <- apply(X,1,sum)
  this_U = rep(1, nrow(X)) #initialize U as a constant vector with mean 1
  this_V = colSums(X)/colSums(L) #initialize V as the mean of the data

  U_record = matrix(data = 0, nrow = iter, ncol = length(this_U), byrow = TRUE)
  V_record = matrix(data = 0, nrow = iter, ncol = length(this_V), byrow = TRUE)
  #pii_record = rep(0, iter)
  U_record[1, ] = rowSums(X)
  V_record[1, ] = this_V

  for(m in 2:iter){

    ptm <- proc.time()
    ### UtoV
    this_V <- Amat_sumi / this_U %*% L
    ### VtoU
    this_U <- Amat_sumj / this_V %*% t(L)

    U_record[m, ] = this_U
    V_record[m, ] = this_V

    print(c(round(m), (proc.time() - ptm)[3]))
    ptm <- proc.time()
  }
  return(list(U_record, V_record))
}

ptm <- proc.time()
ATAC_vec = vector_generator(X, L, iter)
#DNase_vec = vector_generator_simplex(Y, L, iter)

U <- t(ATAC_vec[[1]])
rownames(U) <- rownames(X)
#B <- t(ATAC_vec[[2]])
V <- t(ATAC_vec[[2]])
#PI <- ATAC_vec[[3]]

write.table(round(U,4),file=paste(tissue,flankN,datatype,"U_iterPos.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(round(V,4),file=paste(tissue,flankN,datatype,"V_iterPos.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(round(B,4),file=paste(tissue,flankN,datatype,"B_simplexZIP.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(round(PI,4),file=paste(tissue,flankN,datatype,"PI_ZIP.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
outU <- cbind(U[,1],U[,ncol(U)])
rownames(outU) <- rownames(U)
write.table(round(outU,4),file=paste(tissue,flankN,datatype,"outU_iterPos.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)


if(as.numeric(flankN) == 1){
  PCH <- 16
}else{PCH <- "."}

corpeak <- cor(U[,1],U[,2:ncol(U)])[1,]

selfcorU <- c()
for(i in 2:(ncol(U)-1)){
  selfcorU <- c(selfcorU,cor(U[,i],U[,i+1]))
}

selfcorV <- c()
for(i in 2:(ncol(V)-1)){
  selfcorV <- c(selfcorV,cor(V[,i],V[,i+1]))
}

#selfcorB <- c()
#for(i in 2:(ncol(B)-1)){
#  selfcorB <- c(selfcorB,cor(B[,i],B[,i+1]))
#}

pdf(file=paste(tissue,flankN,datatype,"_iterPos.pdf",sep="_"),height=9,width=9)
par(mfrow=c(3,3),mar=c(4,4,2,2))
#plot(as.numeric(PI),type="h",xlab="iterN",ylab="Pi(% from Zero)",main=paste(tissue,flankN,datatype,sep=" "))
#axis(side=1)
plot(1,1,type="n",xlab="",ylab="",axes=F,main=paste(tissue,flankN,datatype,sep=" "))

plot(corpeak,type="l",xlab="iterN",ylab="Peaksig cor U",main="Peaksig cor U")
plot(selfcorU,type="l",xlab="iterN",ylab="U self-cor",main="U self-cor")
plot(selfcorV,type="l",xlab="iterN",ylab="V self-cor",main="V self-cor")
#plot(selfcorB,type="l",xlab="iterN",ylab="B self-cor",main="B self-cor")
plot(1,1,type="n",xlab="",ylab="",axes=F,main="")

plot(outU[,1], outU[,2],pch=".",xlab="Peaksig",ylab="lastU",main="Peaksig vs lastU")
legend("topleft",legend=paste0("C.C = ",round(cor(outU[,1], outU[,2]),4)),bty="n")
plot(exp(bias),V[,ncol(V)],pch=PCH,xlab="encbias",ylab="lastV",main="encbias vs lastV")
legend("topleft",legend=paste0("C.C = ",round(cor(exp(bias), V[,ncol(V)]),4)),bty="n")
#plot(simplex[2:length(simplex)], apply(B[,(ncol(B)-100):ncol(B)],1,mean),pch=16,xlab="encParam",ylab="last100B",main="encbias vs last100B")
#legend("topleft",legend=paste0("C.C = ",round(cor(simplex[2:length(simplex)], apply(B[,(ncol(B)-100):ncol(B)],1,mean)),4)),bty="n")
plot(1,1,type="n",xlab="",ylab="",axes=F,main="")

ATAC_ratio <- as.numeric(outU[,2]/(outU[,1]+1))
names(ATAC_ratio)<-rownames(outU)
boxplot(ATAC_ratio,ATAC_ratio[TSSpeak],ATAC_ratio[TTSpeak],ATAC_ratio[otherpeak],outline=F,names=c("all","TSS","TTS","other"),ylab="U/reads",las=2)
ATTS <- -log10(t.test(ATAC_ratio[TTSpeak],ATAC_ratio,alternative="less")$p.val)
ATSS <- -log10(t.test(ATAC_ratio[TSSpeak],ATAC_ratio,alternative="less")$p.val)
Aother <- -log10(t.test(ATAC_ratio[otherpeak],ATAC_ratio,alternative="less")$p.val)
legend("topleft",legend=paste0(c("TSS","TTS","other"),":",c(ATSS,ATTS,Aother)),bty="n")
#return(c(ATSS,ATTS,Aother))
dev.off()

