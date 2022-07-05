a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- as.numeric(a[2])#1
datatype <- a[3]#"ATAC"
iter = 1000

ATAC_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_",datatype,"_readscount_seqtype.bed"),header=T,row.names=4)
seqtype_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T,row.names=4)
Hmat <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Data/bias_matrix/simplex_Hmat/simplex_",flankN*2,"mer_Hmat.txt"),row.names=1)

GC_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_ADmergePeak_summit200_GC.bed"),row.names=4,header=F)

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
G <- GC_raw[,4]

X = matrix(as.numeric(as.matrix(X)), ncol = ncol(X))
L = matrix(as.numeric(as.matrix(L)), ncol = ncol(L))
G = as.numeric(G)/400

rownames(X) <- rownames(ATAC_raw)
rownames(L) <- rownames(ATAC_raw)
names(G) <- rownames(ATAC_raw)

pos = which(rowSums(X) > 0 & rowSums(X) < 1000)

X = X[pos, ]
L = L[pos, ]
G = G[pos]
H <- matrix(as.numeric(as.matrix(Hmat)), ncol = ncol(Hmat))[,2:ncol(Hmat)]


#G = runif(nrow(X), 0.5, 1.5)

vector_generator_GC = function(X, L, G, iter){
  #in this version, lambda = (exp(U) + G) * exp(H %*% B) * L
  set.seed(1228)
  #initial the parameters
  likeli = rep(0, iter) #likelihood
  I = rep(0, ncol(X) * nrow(X)) #indicator to show whether it is from Poisson distribution
  I[which(X > 0)] = 1
  I[which(L == 0)] = 0
  pos = intersect(which(X == 0), which(L > 0))
  I[pos] = rbinom(length(pos), 1, 0.5) #initial I is random drawn
  I = matrix(I, ncol = ncol(X))
  #U = runif(nrow(X)) #initial U is from bounded uniform(0,1)
  U = runif(nrow(X), -0.5, 0.5) #initial U are random drawn from (-0.5, 0.5)
  B = runif(ncol(H), -0.5, 0.5) #initial B are random drawn from (-0.5, 0.5)
  # V = exp(H %*% B)
  pii = rbeta(1, 1 + length(which(L > 0)) - sum(I), 1 + sum(I)) #initial pii is from beta distribution
  
  #matrix/vector to record U, V and pi
  U_record = matrix(data = 0, nrow = iter, ncol = length(U), byrow = TRUE)
  W_record = matrix(data = 0, nrow = iter, ncol = length(U), byrow = TRUE)
  B_record = matrix(data = 0, nrow = iter, ncol = length(B), byrow = TRUE)
  V_record = matrix(data = 0, nrow = iter, ncol = nrow(H), byrow = TRUE)
  pii_record = rep(0, iter)
  U_record[1, ] = U
  W_record[1, ] = exp(U)
  B_record[1, ] = B
  V_record[1, ] = exp(H %*% B)
  pii_record[1] = pii
  
  #start the iteration
  for(m in 2:iter){
    V = exp(H %*% B)
    # IX = I * X
    # IL = I * L
    #update U
    U_new = U
    temp_U = runif(length(U), -0.1, 0.1) #random walk of U
    U_new = U + temp_U
    # pos_rej = union(which(U_new <= 0), which(U_new >= 2))
    # U_new[pos_rej] = U[pos_rej] #reject those out of range
    #calculate the difference of 
    diff = rowSums((X == 0) * (L > 0) * (log(pii + (1 - pii) * exp(-(exp(U_new) + G) %*% t(V) * L)) - 
                                 log(pii + (1 - pii) * exp(-(exp(U) + G) %*% t(V) * L)))) +
      rowSums(X) * (log(exp(U_new) + G) - log(exp(U) + G)) - rowSums((X != 0) * ((exp(U_new) - exp(U)) %*% t(V) * L))

    pos_change = which(runif(length(U)) < exp(diff))
    U[pos_change] = U_new[pos_change]
    
    #update B
    B_new = B
    temp_B = runif(length(B), -0.1, 0.1) #random walk of B
    B_new = B + temp_B
    # pos_rej = union(which(B_new <= -10), which(B_new >= 10))
    # B_new[pos_rej] = B[pos_rej]
    # temp_B[pos_rej] = 0 #reject those out of range
    mat_B = t(matrix(rep(temp_B, nrow(H)), nrow = length(temp_B)))
    mod_U = exp(U) + G
    V = exp(H %*% B)
    diff = colSums(log(pii + (1 - pii) * exp((-L * (X == 0) * (L > 0) * mod_U %*% t(V)) %*% exp(mat_B * H)))) - 
      rep(sum((X == 0) * (L > 0) * log(pii + (1 - pii) * exp(-L * mod_U %*% t(V)))), length(temp_B)) + 
      colSums(X %*% H) * temp_B - colSums(((X != 0) * L * (mod_U %*% t(V))) %*% (exp(mat_B * H) - 1))
    
    # diff = temp_B * colSums(IX %*% H) - colSums((IL * (U %*% t(V))) %*% (exp(mat_B * H) - 1))
    pos_change = which(runif(length(B)) < exp(diff))
    B[pos_change] = B_new[pos_change]
    
    #record U, B and V
    U_record[m, ] = U
    W_record[m, ] = exp(U)
    B_record[m, ] = B
    V_record[m, ] = exp(H %*% B)
    
    #update I
    J = as.vector(pii/(pii + (1 - pii) * exp(-(exp(U) + G) %*% t(exp(H %*% B)) * L)))
    I = runif(length(J)) > J
    I[which(X > 0)] = 1
    I[which(L == 0)] = 0
    I = matrix(I, ncol = ncol(X))
    
    #update pii
    pii = rbeta(1, 1 + length(which(L > 0)) - sum(I), 1 + sum(I))
    pii_record[m] = pii
    print(m)
    
    lambda_esi = (exp(U) + G) %*% t(exp(H %*% B)) * L
    pos_zero = intersect(which(X == 0), which(L > 0))
    pos_non = which(X != 0)
    l_zero = log(pii + (1 - pii) * exp(-lambda_esi[pos_zero]))
    l_non = -lambda_esi[pos_non] + X[pos_non] * log(lambda_esi[pos_non])
    likeli[m] = sum(l_zero) + sum(l_non)
    print(likeli[m])
    print(c(round(m), round(pii,3), round(likeli[m],3), (proc.time() - ptm)[3]))
    ptm <- proc.time()
  }
  return(list(U_record, W_record, B_record, V_record, pii_record, likeli))
}

ptm <- proc.time()
ATAC_vec = vector_generator_GC(X, L, G, iter)

U <- t(ATAC_vec[[1]])
rownames(U) <- rownames(X)
W <- t(ATAC_vec[[2]])
B <- t(ATAC_vec[[3]])
V <- t(ATAC_vec[[4]])
PI <- ATAC_vec[[5]]
LIKELI <- ATAC_vec[[6]]


write.table(round(U,4),file=paste(tissue,flankN,datatype,"U_simplexZIPlikeliGC.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(round(W,4),file=paste(tissue,flankN,datatype,"W_simplexZIPlikeliGC.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(round(V,4),file=paste(tissue,flankN,datatype,"V_simplexZIPlikeliGC.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(round(B,4),file=paste(tissue,flankN,datatype,"B_simplexZIPlikeliGC.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(round(PI,4),file=paste(tissue,flankN,datatype,"PI_simplexZIPlikeliGC.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(round(LIKELI,4),file=paste(tissue,flankN,datatype,"likeli_simplexZIPlikeliGC.txt",sep="_"),row.names=F,col.names=F,sep="\t",quote=F)
outU <- cbind(U[,1],apply(U[,(ncol(U)-100):ncol(U)],1,mean))
rownames(outU) <- rownames(U)
write.table(round(outU,4),file=paste(tissue,flankN,datatype,"outU_simplexZIPlikeliGC.txt",sep="_"),row.names=T,col.names=F,sep="\t",quote=F)



#ATAC_vec = vector_generator_GC(X, L, H, G, iter)
#
#DNase_vec = vector_generator_GC(Y, L, H, G, iter)
#
#uni = sort(unique(rowSums(X)))
#sdd = rep(0, length(uni))
#for(i in 1:length(uni)){
#  temp = which(rowSums(X) == uni[i])
#  sdd[i] = sd(ATAC_vec[[1]][500, temp])/mean(ATAC_vec[[1]][500, temp])
#}
#
#plot(rowSums(X), colSums(ATAC_vec[[1]][300:500, ]))
#plot(rowSums(Y), colSums(DNase_vec[[1]][300:500, ]))
#plot(rowSums(Y), colSums(DNase_vec[[1]][450:500, ]))
#
#
#cor(rowSums(X), rowSums(Y)) #original correlation
#cor(colSums(ATAC_vec[[1]]), rowSums(X))
#cor(colSums(DNase_vec[[1]]), rowSums(Y))
#cor(colSums(ATAC_vec[[1]]), colSums(DNase_vec[[1]])) #correlation of generated vector, history mean
#cor(colSums(ATAC_vec[[1]][300:500,]), colSums(DNase_vec[[1]][300:500,])) #correlation of generated vector, history mean
#cor(colSums(ATAC_vec[[1]][300:500,]), rowSums(X))
#cor(colSums(DNase_vec[[1]][300:500,]), rowSums(Y))
#
#write.table(ATAC_vec[[1]], file = "U_ATAC4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(ATAC_vec[[2]], file = "B_ATAC4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(ATAC_vec[[3]], file = "V_ATAC4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(ATAC_vec[[4]], file = "pi_ATAC4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(DNase_vec[[1]], file = "U_DNase4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(DNase_vec[[2]], file = "H_DNase4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(DNase_vec[[3]], file = "V_DNase4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(DNase_vec[[4]], file = "pi_DNase4_Ovary.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#
#
#
#
#