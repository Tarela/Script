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
  set.seed(1228)

  #initialize the iteration
  #I = rep(0, ncol(X) * nrow(X))
  #I[which(X > 0)] = 1
  #I[which(L == 0)] = 0
  pos = intersect(which(X == 0), which(L > 0))
  #I[pos] = rbinom(length(pos), 1, 0.5) #initial I is random drawn
  #I = matrix(I, ncol = ncol(X))
  # U = runif(nrow(X)) #initial U is a random drawn
  # V = runif(ncol(X)) #initial V is a random drawn
  U = rep(1, nrow(X)) #initialize U as a constant vector with mean 1
  V = colSums(X)/colSums(L) #initialize V as the mean of the data
  #pii = rbeta(1, 1+length(which(L>0))-sum(I), 1+sum(I)) #initial pii is from beta distribution
  
  #matrix/vector to record U, V and pi
  U_record = matrix(data = 0, nrow = iter, ncol = length(U), byrow = TRUE)
  V_record = matrix(data = 0, nrow = iter, ncol = length(V), byrow = TRUE)
  #pii_record = rep(0, iter)
  U_record[1, ] = rowSums(X)
  V_record[1, ] = V
  #pii_record[1] = pii
  #pi_record[1] = 1 - sum(I)/length(which(L > 0))
  temp_U = rep(0, length(U))
  temp_V = rep(0, length(V))
  #since we have to give a proper prior, we need to find a bound for U * exp(B %*% H)
  pos_L = L[which(L > 0)]
  pos_X = X[which(L > 0)]
  bound = max(pos_X/pos_L)
  
  #then begin the iteration
  for(k in 2:iter){
    #update U
    for(i in 1:length(U)){
      U_new = U
      temp_U = runif(1, -0.05, 0.05) #random walk of U
      U_new[i] = U[i] + temp_U #sequentailly update U
#      if(U_new[i] >= 0 && U_new[i] <= sqrt(bound)){
      if(U_new[i] >= 0 && U_new[i] <= 2){
        diff = sum(X[i, ] * (log(U_new[i]) - log(U[i]))) - temp_U * sum(V * L[i, ])
        if(runif(1) < exp(diff)){
          U[i] = U_new[i] #U is updated
        }
      }
    }
    #update V
    for(j in 1:length(V)){
      V_new = V
      temp_V = runif(1, -0.05, 0.05) 
      V_new[j] = V[j] + temp_V #random walk of V bounded in [0, sqrt(bound)]
      if(V_new[j] >= 0 && V_new[j] <= sqrt(bound)){
        diff = sum(X[, j] * (log(V_new[j]) - log(V[j]))) - temp_V * sum(U * L[, j])
        if(runif(1) < exp(diff)){
          V[j] = V_new[j] #V is updated
        }
      }
    }
    # #normalize U and V. The mean of U remains 1 and keep the multiply of U and V unchanged
    # V = V * mean(U)
    # U = U / mean(U)
    U_record[k, ] = U
    V_record[k, ] = V
    #update I
    #J = as.vector(pii/(pii + (1 - pii) * exp(-U %*% t(V) * L)))
    #I = runif(length(J)) > J
    #I[which(X > 0)] = 1
    #I[which(L == 0)] = 0
    #I = matrix(I, ncol = ncol(X))
    #update pii
    #pii = rbeta(1, 1+length(which(L>0))-sum(I), 1+sum(I))
    #pii_record[k] = pii
    print(c(round(k), (proc.time() - ptm)[3]))
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

write.table(round(U,4),file=paste(tissue,flankN,datatype,"U_samplingPos.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(round(V,4),file=paste(tissue,flankN,datatype,"V_samplingPos.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(round(B,4),file=paste(tissue,flankN,datatype,"B_simplexZIP.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(round(PI,4),file=paste(tissue,flankN,datatype,"PI_ZIP.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
outU <- cbind(U[,1],apply(U[,(ncol(U)-100):ncol(U)],1,mean))
rownames(outU) <- rownames(U)
write.table(round(outU,4),file=paste(tissue,flankN,datatype,"outU_samplingPos.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)


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

pdf(file=paste(tissue,flankN,datatype,"_samplingPos.pdf",sep="_"),height=9,width=9)
par(mfrow=c(3,3),mar=c(4,4,2,2))
#plot(as.numeric(PI),type="h",xlab="iterN",ylab="Pi(% from Zero)",main=paste(tissue,flankN,datatype,sep=" "))
#axis(side=1)
plot(1,1,type="n",xlab="",ylab="",axes=F,main=paste(tissue,flankN,datatype,sep=" "))

plot(corpeak,type="l",xlab="iterN",ylab="Peaksig cor U",main="Peaksig cor U")
plot(selfcorU,type="l",xlab="iterN",ylab="U self-cor",main="U self-cor")
plot(selfcorV,type="l",xlab="iterN",ylab="V self-cor",main="V self-cor")
#plot(selfcorB,type="l",xlab="iterN",ylab="B self-cor",main="B self-cor")
plot(1,1,type="n",xlab="",ylab="",axes=F,main="")

plot(outU[,1], outU[,2],pch=".",xlab="Peaksig",ylab="last100U",main="Peaksig vs last100U")
legend("topleft",legend=paste0("C.C = ",round(cor(outU[,1], outU[,2]),4)),bty="n")
plot(exp(bias), apply(V[,(ncol(V)-100):ncol(V)],1,mean),pch=PCH,xlab="encbias",ylab="last100V",main="encbias vs last100V")
legend("topleft",legend=paste0("C.C = ",round(cor(exp(bias), apply(V[,(ncol(V)-100):ncol(V)],1,mean)),4)),bty="n")
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
