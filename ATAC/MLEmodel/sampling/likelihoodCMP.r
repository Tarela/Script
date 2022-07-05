a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- as.numeric(a[2])#1
datatype <- a[3]#"ATAC"
#iter = 1000


U_sampling_raw <- read.table(paste0(paste(tissue,flankN,datatype,"outU","samplingPos",sep="_"),".txt"), row.names=1)
U_GD_raw <- read.table(paste0(paste(tissue,flankN,datatype,"outU","GDsimplex",sep="_"),".txt"), row.names=1)
V_sampling_raw <- read.table(paste0(paste(tissue,flankN,datatype,"V","samplingPos",sep="_"),".txt"))
V_GD_raw <- read.table(paste0(paste(tissue,flankN,datatype,"V","GDsimplex",sep="_"),".txt"))

U_sampling <-  U_sampling_raw[,2]
U_GD <-  U_GD_raw[,2]
V_sampling <-  apply(V_sampling_raw[,900:1000],1,mean)
V_GD <-  V_GD_raw[,ncol(V_GD_raw)]

ATAC_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_",datatype,"_readscount_seqtype.bed"),header=T,row.names=4)
seqtype_raw <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/simplexZIP_application/CAR_features/peak_matrix/ENCODEtissue/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T,row.names=4)


X <- ATAC_raw[,4:ncol(ATAC_raw)]#[1:1000,]
L <- seqtype_raw[,4:ncol(seqtype_raw)]#[1:1000,]

X = matrix(as.numeric(as.matrix(X)), ncol = ncol(X))
L = matrix(as.numeric(as.matrix(L)), ncol = ncol(L))
rownames(X) <- rownames(ATAC_raw)
rownames(L) <- rownames(ATAC_raw)
#pos = which(rowSums(X) > 0 & rowSums(X) < 1000)

X = X[rownames(U_sampling_raw), ]
L = L[rownames(U_sampling_raw), ]
#H <- matrix(as.numeric(as.matrix(Hmat)), ncol = ncol(Hmat))[,2:ncol(Hmat)]

#useU <- U[,ncol(U)]
#useV <- V[,ncol(V)]
#
sumlikeli <- function(U,V,X,L){
  UVL <- as.numeric((U %*% t(V))*L)
  useX <- as.numeric(X)
  idx <- which(UVL>0)
  likeli_single <- useX[idx]*log(UVL[idx]) - UVL[idx]
  return(sum(likeli_single))
}

print(c(tissue,flankN,datatype,"sampling",sumlikeli(U_sampling,V_sampling,X,L)))
print(c(tissue,flankN,datatype,"GD",sumlikeli(U_GD,V_GD,X,L)))
