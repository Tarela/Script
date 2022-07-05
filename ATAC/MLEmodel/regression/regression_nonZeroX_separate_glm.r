#a <- c(1:100)
#y <- rnorm(100)
#
#lm(y ~ a - 1)
#lm(y ~ a)
#lm(y ~ as.factor(a))
#lm(y ~ as.factor(a) - 1)
#



### prepare data:
a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- a[2]#1
datatype <- a[3]#"ATAC"

ATAC_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount_new/",tissue,"_mergePeak_flank",flankN,"_",datatype,"_readscount_seqtype.bed"),header=T)
seqtype_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T)


ATAC_mat <- ATAC_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(ATAC_raw)]#[1:100,]
seqtype_mat <- seqtype_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(seqtype_raw)]#[1:100,]

## dim
M <- nrow(ATAC_mat)
N <- ncol(ATAC_mat)

X <- rapply(ATAC_mat,c)
L <- rapply(seqtype_mat,c)

alpha <- rep(seq(M),N)
beta <- rep(seq(N),each=M)


idx <- which(X>0)
use_X <- X[idx]
use_L <- L[idx]
use_alpha <- alpha[idx]
use_beta <- beta[idx]

#reg_model <- glm(X~as.factor(alpha)+as.factor(beta),family="poisson")
reg_model <- glm(use_X~as.factor(use_alpha)+as.factor(use_beta)+use_L,family="poisson")
#reg_model <- lm(X~as.factor(alpha)+as.factor(beta))

reg_coeff <- reg_model$coefficients
### 
Intercept <- reg_coeff["(Intercept)"]
coeffL <- reg_coeff["use_L"]
coeffA <- c(0,reg_coeff[2:M])
coeffB <- c(0,reg_coeff[(M+1):(M+N-1)])

logV <- rep(0,N)
logV[1] <-  -1*sum(coeffB,na.rm=T)/N 
for(j in 2:N){
    coeffB_name <- paste0("as.factor(use_beta)",j)
    if (coeffB_name %in% names(coeffB)){
        logV[j] <- coeffB[coeffB_name]
    }
#    logV <- c(logV,coeffB[j]+logV[1])
}

logU <- rep(0,M)
logU[1] <- Intercept + coeffL*L[1] - logV[1]
for(i in 2:M){
    coeffA_name <- paste0("as.factor(use_alpha)",i)
    if (coeffA_name %in% names(coeffA)){
        logU[i] <- coeffA[coeffA_name]
    }
    #logU <- c(logU,coeffA[i]+logU[1])
}

ATAC_counts <- apply(ATAC_mat,1,sum)

#colnames(X)
names(logV) <- colnames(ATAC_mat)

write.table(x=logV,file=paste0(tissue,"_mergePeak_flank",flankN,"_",datatype,"_biasV.txt"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(x=cbind(ATAC_counts, logU),file=paste0(tissue,"_mergePeak_flank",flankN,"_",datatype,"_countU.txt"),row.names=F,col.names=F,sep="\t",quote=F)

#test <- cbind(X,alpha,beta,L)
#U <- matrix(rep(0,M*N*M), nrow=M*N)





