#a <- c(1:100)
#y <- rnorm(100)
#
#lm(y ~ a - 1)
#lm(y ~ a)
#lm(y ~ as.factor(a))
#lm(y ~ as.factor(a) - 1)
#

library(biglm)


### prepare data:
a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- a[2]#1
#datatype <- a[3]#"ATAC"

ATAC_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_ATAC_readscount_seqtype.bed"),header=T)
DNase_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_DNase_readscount_seqtype.bed"),header=T)
seqtype_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T)


ATAC_mat <- ATAC_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(ATAC_raw)]#[1:100,]
DNase_mat <- DNase_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(ATAC_raw)]#[1:100,]
seqtype_mat <- seqtype_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(seqtype_raw)]#[1:100,]

## dim
#combine_mat <- rbind(ATAC_mat,DNase_mat)
#combine_seqtype_mat <- rbind(seqtype_mat,seqtype_mat)

M <- nrow(ATAC_mat)
N <- ncol(ATAC_mat)

#X <- rapply(combine_mat,c)
XA <- rapply(ATAC_mat,c)
XD <- rapply(DNase_mat,c)
X <- c(XA,XD)
L <- c(rapply(seqtype_mat,c),rapply(seqtype_mat,c))

alpha <- rep(seq(M),2*N)
beta <- rep(seq(N*2),each=M)
#betaD <- rep(seq(N),each=M)
#reg_model <- glm(X~as.factor(alpha)+as.factor(beta),family="poisson")
#reg_model <- glm(X~as.factor(alpha)+as.factor(beta)+L,family="poisson")
#reg_model <- lm(X~as.factor(alpha)+as.factor(beta))
#reg_coeff <- reg_model$coefficients

usedata <- as.data.frame(cbind(X,alpha,beta,L))
usedata$alpha_factor <- factor(usedata$alpha)
usedata$beta_factor <- factor(usedata$beta)
#reg_model1 <- glm(X~as.factor(alpha)+as.factor(beta)+L,data=usedata,family=poisson(log))
#reg_model <- lm(X~as.factor(alpha)+as.factor(beta))
reg_model <- bigglm(X~alpha_factor+beta_factor+L,data=usedata,family=poisson(log),maxit=100,chunksize=10000)
#reg_model <- bigglm(X~as.factor(alpha)+as.factor(beta)+L,data=usedata,family=poisson(log),maxit=100)
#reg_model1 <- glm(X~as.factor(alpha)+as.factor(beta)+L,data=usedata,family=poisson(log),maxit=100)

reg_coeff <- summary(reg_model)$mat[,"Coef"]

### 
Intercept <- reg_coeff["(Intercept)"]
coeffL <- reg_coeff["L"]
coeffA <- c(0,reg_coeff[2:M])
coeffB <- c(0,reg_coeff[(M+1):(M+2*N-1)])

logV <- c( -1*sum(coeffB)/N )
for(j in 2:(2*N)){
    logV <- c(logV,coeffB[j]+logV[1])
}

logU <- c( Intercept + coeffL*L[1] - logV[1])
for(i in 2:M){
    logU <- c(logU,coeffA[i]+logU[1])
}

ATAC_counts <- apply(ATAC_mat,1,sum)
DNase_counts <- apply(DNase_mat,1,sum)
logV_ATAC <- logV[1:(4**(2*flankN))]
logV_DNase <- logV[(4**(2*flankN)+1):((4**(2*flankN))*2)]
#colnames(X)
names(logV_ATAC) <- colnames(ATAC_mat)
names(logV_DNase) <- colnames(ATAC_mat)

write.table(x=cbind(logV_ATAC,logV_DNase),file=paste0(tissue,"_mergePeak_flank",flankN,"_combine_biasV.txt"),row.names=T,col.names=F,sep="\t",quote=F)
write.table(x=cbind(ATAC_counts, DNase_counts, logU),file=paste0(tissue,"_mergePeak_flank",flankN,"_combine_countU.txt"),row.names=F,col.names=F,sep="\t",quote=F)

#test <- cbind(X,alpha,beta,L)
#U <- matrix(rep(0,M*N*M), nrow=M*N)





