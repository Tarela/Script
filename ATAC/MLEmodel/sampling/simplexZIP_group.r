a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- as.numeric(a[2])#1
datatype <- a[3]#"ATAC"
iter = 1000

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


vector_generator_simplex = function(X, L, iter){
  set.seed(1228)

  I = rep(0, ncol(X) * nrow(X))
  I[which(X > 0)] = 1
  I[which(L == 0)] = 0
  pos = intersect(which(X == 0), which(L > 0))
  I[pos] = rbinom(length(pos), 1, 0.5) #initial I is random drawn
  I = matrix(I, ncol = ncol(X))
  #U = runif(nrow(X)) #initial U is from bounded uniform(0,1)
  U = rep(1, nrow(X)) #initial U are all 1, and bound it into [0,2]
  B = runif(ncol(H), -1, 1) #initial B is random drawn from (-1, 1)
  # V = exp(H %*% B)
  pii = rbeta(1, 1+length(which(L>0))-sum(I), 1+sum(I)) #initial PI is from beta distribution
   
  #matrix/vector to record U, V and pi
  U_record = matrix(data = 0, nrow = iter, ncol = length(U), byrow = TRUE)
  B_record = matrix(data = 0, nrow = iter, ncol = length(B), byrow = TRUE)
  V_record = matrix(data = 0, nrow = iter, ncol = nrow(H), byrow = TRUE)
  pii_record = rep(0, iter)
  U_record[1, ] = rowSums(X)
  B_record[1, ] = B
  V_record[1, ] = exp(H %*% B)
  pii_record[1] = pii
  
  #start the iteration
  for(m in 2:iter){
    V = exp(H %*% B)
    IX = I * X
    IL = I * L
    #update U
    U_new = U
    temp_U = runif(length(U), -0.05, 0.05) #random walk of U
    U_new = U + temp_U
    pos_rej = union(which(U_new <= 0), which(U_new >= 2))
    U_new[pos_rej] = U[pos_rej] #reject those out of range
    diff = rowSums(IX) * (log(U_new) - log(U)) - (U_new - U) * (IL %*% V)
    pos_change = which(runif(length(U)) < exp(diff))
    U[pos_change] = U_new[pos_change]

    #update B
    B_new = B
    temp_B = runif(length(B), -0.05, 0.05) #random walk of B
    B_new = B + temp_B
    pos_rej = union(which(B_new <= -10), which(B_new >= 10))
    B_new[pos_rej] = B[pos_rej]
    temp_B[pos_rej] = 0 #reject those out of range
    mat_B = t(matrix(rep(temp_B, nrow(H)), nrow = length(temp_B)))
    diff = temp_B * colSums(IX %*% H) - colSums((IL * (U %*% t(V))) %*% (exp(mat_B * H) - 1))
    pos_change = which(runif(length(B)) < exp(diff))
    B[pos_change] = B_new[pos_change]
    
    #record U, B and V
    U_record[m, ] = U
    B_record[m, ] = B
    V_record[m, ] = exp(H %*% B)
    
    #update I
    J = as.vector(pii/(pii + (1 - pii) * exp(-U %*% t(exp(H %*% B)) * L)))
    I = runif(length(J)) > J
    I[which(X > 0)] = 1
    I[which(L == 0)] = 0
    I = matrix(I, ncol = ncol(X))
    
    #update pii
    pii = rbeta(1, 1+length(which(L>0))-sum(I), 1+sum(I))
    pii_record[m] = pii
    print(c(round(m), round(pii,3), (proc.time() - ptm)[3]))
    ptm <- proc.time()
  }
  return(list(U_record, B_record, V_record, pii_record))
}

ptm <- proc.time()
ATAC_vec = vector_generator_simplex(X, L, iter)
#DNase_vec = vector_generator_simplex(Y, L, iter)

U <- t(ATAC_vec[[1]])
rownames(U) <- rownames(X)
B <- t(ATAC_vec[[2]])
V <- t(ATAC_vec[[3]])
PI <- ATAC_vec[[4]]

write.table(round(U,4),file=paste(tissue,flankN,datatype,"U_simplexZIPgroup.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(round(V,4),file=paste(tissue,flankN,datatype,"V_simplexZIPgroup.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(round(B,4),file=paste(tissue,flankN,datatype,"B_simplexZIPgroup.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(round(PI,4),file=paste(tissue,flankN,datatype,"PI_simplexZIPgroup.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
outU <- cbind(U[,1],apply(U[,(ncol(U)-100):ncol(U)],1,mean))
rownames(outU) <- rownames(U)
write.table(round(outU,4),file=paste(tissue,flankN,datatype,"outU_simplexZIPgroup.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)


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

selfcorB <- c()
for(i in 2:(ncol(B)-1)){
  selfcorB <- c(selfcorB,cor(B[,i],B[,i+1]))
}

pdf(file=paste(tissue,flankN,datatype,"_simplexZIPgroup.pdf",sep="_"),height=9,width=9)
par(mfrow=c(3,3),mar=c(4,4,2,2))
plot(as.numeric(PI),type="h",xlab="iterN",ylab="Pi(% from Zero)",main=paste(tissue,flankN,datatype,sep=" "))
axis(side=1)

plot(corpeak,type="l",xlab="iterN",ylab="Peaksig cor U",main="Peaksig cor U")
plot(selfcorU,type="l",xlab="iterN",ylab="U self-cor",main="U self-cor")
plot(selfcorV,type="l",xlab="iterN",ylab="V self-cor",main="V self-cor")
plot(selfcorB,type="l",xlab="iterN",ylab="B self-cor",main="B self-cor")

plot(outU[,1], outU[,2],pch=".",xlab="Peaksig",ylab="last100U",main="Peaksig vs last100U")
legend("topleft",legend=paste0("C.C = ",round(cor(U[,1], apply(U[,(ncol(U)-100):ncol(U)],1,mean)),4)),bty="n")
plot(exp(bias), apply(V[,(ncol(V)-100):ncol(V)],1,mean),pch=PCH,xlab="encbias",ylab="last100V",main="encbias vs last100V")
legend("topleft",legend=paste0("C.C = ",round(cor(exp(bias), apply(V[,(ncol(V)-100):ncol(V)],1,mean)),4)),bty="n")
plot(simplex[2:length(simplex)], apply(B[,(ncol(B)-100):ncol(B)],1,mean),pch=16,xlab="encParam",ylab="last100B",main="encbias vs last100B")
legend("topleft",legend=paste0("C.C = ",round(cor(simplex[2:length(simplex)], apply(B[,(ncol(B)-100):ncol(B)],1,mean)),4)),bty="n")

ATAC_ratio <- as.numeric(outU[,2]/(outU[,1]+1))
names(ATAC_ratio)<-rownames(outU)
boxplot(ATAC_ratio,ATAC_ratio[TSSpeak],ATAC_ratio[TTSpeak],ATAC_ratio[otherpeak],outline=F,names=c("all","TSS","TTS","other"),ylab="U/reads",las=2)
ATTS <- -log10(t.test(ATAC_ratio[TTSpeak],ATAC_ratio,alternative="less")$p.val)
ATSS <- -log10(t.test(ATAC_ratio[TSSpeak],ATAC_ratio,alternative="less")$p.val)
Aother <- -log10(t.test(ATAC_ratio[otherpeak],ATAC_ratio,alternative="less")$p.val)
legend("topleft",legend=paste0(c("TSS","TTS","other"),":",c(ATSS,ATTS,Aother)),bty="n")
#return(c(ATSS,ATTS,Aother))
dev.off()








