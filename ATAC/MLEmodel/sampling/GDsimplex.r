a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- as.numeric(a[2])#1
datatype <- a[3]#"ATAC"
#iter = 100

ATAC_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_",datatype,"_readscount_seqtype.bed"),header=T,row.names=4)
seqtype_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T,row.names=4)
Hmat <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Data/bias_matrix/simplex_Hmat/simplex_",flankN*2,"mer_Hmat.txt"),row.names=1)


if(datatype == "ATAC"){
    bias <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/prediction_separate_reads_cutseq/MLE_model/sampling/biasMat/YeastNaked_ATACPE_enc",flankN,"flank.txt"),row.names=1,header=T)[,2]
    simplex <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Data/bias_matrix/simplex_parameters/ATAC_flank",flankN,"_simplexParam.txt"))[,1]
}else{
    bias <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/prediction_separate_reads_cutseq/MLE_model/sampling/biasMat/IMR90naked_DNaseSE_enc",flankN,"flank.txt"),row.names=1,header=T)[,2]
    simplex <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Data/bias_matrix/simplex_parameters/DNase_flank",flankN,"_simplexParam.txt"))[,1]    
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
H <- matrix(as.numeric(as.matrix(Hmat)), ncol = ncol(Hmat))[,2:ncol(Hmat)]

#initial the hyperparameters, alpha is the step size, and thres is the stopping criterion
alpha = 1e-6
thres = 1
shrink = 0.8
sstep = 0.5

#gradient_method is to use gradient descent method and line search technique to solve the problem
gradient_method = function(X, L, H, alpha, thres, shrink, sstep){
  iter = 1
  dif1 = 2 * thres
  #initialize parameters
  #U = jitter(rowSums(X)/1000)
  U = rep(1, nrow(X))
  B = rep(0, ncol(H))

  U_record = matrix(data = 0, nrow = 1000, ncol = length(U), byrow = TRUE)
  B_record = matrix(data = 0, nrow = 1000, ncol = length(B), byrow = TRUE)
  V_record = matrix(data = 0, nrow = 1000, ncol = nrow(H), byrow = TRUE)
  U_record[iter, ] = rowSums(X)
  B_record[iter, ] = B
  V_record[1, ] = exp(H %*% B)

  temp_U = rep(0, nrow(X))
  temp_B = rep(0, ncol(H))
  #to speed up
  rX = rowSums(X)
  rXH = rowSums(X %*% H)
  while(dif1 > thres){
    likeli = sum(X %*% (H %*% B)) + sum(t(log(U)) %*% X) - sum(U * (L %*% (exp(H %*% B))))
    #temp_U and temp_B are used to find the gredient descent direction
    temp_U = rX/U - L %*% exp(H %*% B)
    for(i in 1:length(B)){
      temp_B[i] = rXH[i] - sum(U * (L %*% (H[,i] * exp(H %*% B))))
    }
    U_new = U + alpha * temp_U
    B_new = B + alpha * temp_B
    likeli_new = sum(X %*% (H %*% B_new)) + sum(t(log(U_new)) %*% X) - sum(U_new * (L %*% (exp(H %*% B_new))))
    dif2 = likeli_new - likeli #if dif2 is large enough, then enter the next iteration, else into the next few lines
    while(dif2 < alpha * sstep * (sum(temp_B * temp_B) + sum(temp_U * temp_U))){
      alpha = alpha * shrink
      U_new = U + alpha * temp_U
      B_new = B + alpha * temp_B
      likeli_new = sum(X %*% (H %*% B_new)) + sum(t(log(U_new)) %*% X) - sum(U_new * (L %*% (exp(H %*% B_new))))
      dif2 = likeli_new - likeli
      #print(c(dif2, alpha, alpha * sstep * (sum(temp_B * temp_B) + sum(temp_U * temp_U))))
    }
    U = U + alpha * temp_U
    B = B + alpha * temp_B
    likeli_new = sum(X %*% (H %*% B)) + sum(t(log(U)) %*% X) - sum(U * (L %*% (exp(H %*% B))))
    dif1 = likeli_new - likeli
    
    iter = iter + 1
    if(iter <= 1000){
      U_record[iter, ] = U
      B_record[iter, ] = B
      V_record[iter, ] = exp(H %*% B)
    }
    #print(paste0("the iteration now is ", iter))
    #print(paste0("the step size is ", alpha))
    #print(paste0("the likelihood in last step is ", likeli))
    #print(paste0("the likelihood in this step is ", likeli_new))
    #print(paste0("the improvement of the likelihood in this step is ", dif1))
    print(c(round(iter), alpha, likeli, likeli_new, dif1, round((proc.time() - ptm)[3],3)))
    ptm <- proc.time()
  }
  return(list(U_record, B_record, V_record, U, B))
}

#use the ATAC and DNase data to get the U vector
ptm <- proc.time()
ATAC_vec = gradient_method(X, L, H, alpha, thres, shrink, sstep)


U <- t(ATAC_vec[[1]])
rownames(U) <- rownames(X)
B <- t(ATAC_vec[[2]])
V <- t(ATAC_vec[[3]])
#lastU <- ATAC_vec[[4]]
#lastB <- ATAC_vec[[5]]
#lastV = exp(H %*% lastB)
idx <- which(apply(U,2,mean)!=0)
U <- U[,idx]
B <- B[,idx]
V <- V[,idx]


write.table(U,file=paste(tissue,flankN,datatype,"U_GDsimplex.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(V,file=paste(tissue,flankN,datatype,"V_GDsimplex.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(B,file=paste(tissue,flankN,datatype,"B_GDsimplex.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(round(PI,4),file=paste(tissue,flankN,datatype,"PI_GDsimplex.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
outU <- cbind(U[,1],U[,ncol(U)])
rownames(outU) <- rownames(U)
write.table(outU,file=paste(tissue,flankN,datatype,"outU_GDsimplex.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)

#useU <- U[,ncol(U)]
#useV <- V[,ncol(V)]
#
#sumlikeli <- function(U,V,X,L){
#  UVL <- as.numeric((U %*% t(V))*L)
#  useX <- as.numeric(X)
#  idx <- which(UVL>0)
#  likeli_single <- useX[idx]*log(UVL[idx]) - UVL[idx]
#  return(sum(likeli_single))
#}
#
#sumlikeli(useU,useV,X,L)







#if(as.numeric(flankN) == 1){
#  PCH <- 16
#}else{PCH <- "."}
#
#corpeak <- c()
#for(i in 2:(ncol(U)-1)){
#  if(mean(U[,i]) != 0){
#    corpeak <- c(corpeak,cor(U[,i],U[,1]))
#  }else{
#    corpeak <- c(corpeak,0)
#  }
#}
#
#selfcorU <- c()
#for(i in 2:(ncol(U)-1)){
#  if(mean(U[,i+1]) != 0){
#    selfcorU <- c(selfcorU,cor(U[,i],U[,i+1]))
#  }
#}
#
#selfcorV <- c()
#for(i in 2:(ncol(V)-1)){
#  if(mean(V[,i+1])!= 0){
#    selfcorV <- c(selfcorV,cor(V[,i],V[,i+1]))
#  }
#}
#
#selfcorB <- c()
#for(i in 2:(ncol(B)-1) ){
#  if(mean(B[,i+1]) != 0){
#    selfcorB <- c(selfcorB,cor(B[,i],B[,i+1]))  
#  }
#}
#
#pdf(file=paste(tissue,flankN,datatype,"_GDsimplex.pdf",sep="_"),height=9,width=9)
#par(mfrow=c(3,3),mar=c(4,4,2,2))
##plot(as.numeric(PI),type="h",xlab="iterN",ylab="Pi(% from Zero)",main=paste(tissue,flankN,datatype,sep=" "))
##axis(side=1)
#plot(1,1,type="n",xlab="",ylab="",axes=F,main=paste(tissue,flankN,datatype,sep=" "))
#
#plot(corpeak,type="l",xlab="iterN",ylab="Peaksig cor U",main="Peaksig cor U")
#plot(selfcorU,type="l",xlab="iterN",ylab="U self-cor",main="U self-cor")
#plot(selfcorV,type="l",xlab="iterN",ylab="V self-cor",main="V self-cor")
#plot(selfcorB,type="l",xlab="iterN",ylab="B self-cor",main="B self-cor")
#
#plot(outU[,1], outU[,2],pch=".",xlab="Peaksig",ylab="lastU",main="Peaksig vs lastU")
#legend("topleft",legend=paste0("C.C = ",round(cor(outU[,1], outU[,2]),4)),bty="n")
#plot(exp(bias), V[,ncol(V)],pch=PCH,xlab="encbias",ylab="lastV",main="encbias vs lastV")
#legend("topleft",legend=paste0("C.C = ",round(cor(exp(bias), V[,ncol(V)]),4)),bty="n")
#plot(simplex[2:length(simplex)],B[,ncol(B)],pch=16,xlab="encParam",ylab="lastB",main="encbias vs lastB")
#legend("topleft",legend=paste0("C.C = ",round(cor(simplex[2:length(simplex)], B[,ncol(B)]),4)),bty="n")
#
#ATAC_ratio <- as.numeric(outU[,2]/(outU[,1]+1))
#names(ATAC_ratio)<-rownames(outU)
#boxplot(ATAC_ratio,ATAC_ratio[TSSpeak],ATAC_ratio[TTSpeak],ATAC_ratio[otherpeak],outline=F,names=c("all","TSS","TTS","other"),ylab="U/reads",las=2)
#ATTS <- -log10(t.test(ATAC_ratio[TTSpeak],ATAC_ratio,alternative="less")$p.val)
#ATSS <- -log10(t.test(ATAC_ratio[TSSpeak],ATAC_ratio,alternative="less")$p.val)
#Aother <- -log10(t.test(ATAC_ratio[otherpeak],ATAC_ratio,alternative="less")$p.val)
#legend("topleft",legend=paste0(c("TSS","TTS","other"),":",c(ATSS,ATTS,Aother)),bty="n")
##return(c(ATSS,ATTS,Aother))
#dev.off()
#
#
#





