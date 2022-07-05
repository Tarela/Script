#a <- c(1:100)
#y <- rnorm(100)
#
#lm(y ~ a - 1)
#lm(y ~ a)
#lm(y ~ as.factor(a))
#lm(y ~ as.factor(a) - 1)
#
## dim
#M <- nrow(ATAC_mat)
#N <- ncol(ATAC_mat)

#X <- rapply(ATAC_mat,c)
#L <- rapply(seqtype_mat,c)

#initU <- (ATAC_counts + DNase_counts)
#initU <- ATAC_counts
#initU <- DNase_counts
#initV <- ATAC_bias[colnames(ATAC_mat),"encbias"]


UtoVW <- function(thisU,Amat,Dmat,Lmat){
    m <- nrow(Amat)
    n <- ncol(Amat)
    newV <- c()
    newW <- c()
    for(j in 1:n){
        Vj = sum(Amat[,j]) / sum(thisU * Lmat[,j])
        newV <- append(newV,Vj)
        Wj = sum(Dmat[,j]) / sum(thisU * Lmat[,j])
        newW <- append(newW,Wj)

    }
    return(cbind(newV,newW))
}

VWtoU <- function(thisVW,Amat,Dmat,Lmat){
    m <- nrow(Amat)
    n <- ncol(Amat)
    #thisV <- thisVW[,1]
    #thisW 
    thisVWcb <- thisVW[,1] + thisVW[,2]
    newU <- c()
    for(i in 1:m){
        Ui = sum(Amat[i,]+Dmat[i,]) / sum(thisVWcb * Lmat[i,])
        newU <- append(newU,Ui)
    }
    return(newU)
}



### prepare data:
a <- commandArgs(T)

tissue <- a[1]#"Ovary"
flankN <- a[2]#1
maxit <- as.numeric(a[3])#100
#datatype <- a[3]#"ATAC"

ATAC_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_ATAC_readscount_seqtype.bed"),header=T)
DNase_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_DNase_readscount_seqtype.bed"),header=T)
seqtype_raw <- read.table(paste0("/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/ATAC_DNase_pairCMP/kmer_readscount_seqcount/",tissue,"_mergePeak_flank",flankN,"_seqcount_seqtype.bed"),header=T)

ATAC_bias <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/diffmer_var_comparison/diffmerMat/flanking/YeastNaked_ATACPE_enc",flankN,"flank.txt"),row.names=1,header=T)#[,"encbias"]
DNase_bias <- read.table(paste0("/scratch/sh8tv/Project/scATAC/Result/diffmer_var_comparison/diffmerMat/flanking/IMR90naked_DNaseSE_enc",flankN,"flank.txt"),row.names=1,header=T)#[,"encbias"]


#ATAC_mat <- ATAC_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(ATAC_raw)]#[1:100,]
#DNase_mat <- DNase_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(ATAC_raw)]#[1:100,]
#seqtype_mat <- seqtype_raw[which(ATAC_raw[,1]=="chr1"),4:ncol(seqtype_raw)]#[1:100,]

ATAC_mat <- ATAC_raw[,4:ncol(ATAC_raw)]#[1:100,]
DNase_mat <- DNase_raw[,4:ncol(ATAC_raw)]
seqtype_mat <- seqtype_raw[,4:ncol(seqtype_raw)]#[1:100,]

ATAC_counts <- apply(ATAC_mat,1,sum)
DNase_counts <- apply(DNase_mat,1,sum)


Amat <- ATAC_mat
Dmat <- DNase_mat
Lmat <- seqtype_mat

summary_U <- cbind(ATAC_counts, DNase_counts)
summary_V <- c(ATAC_bias[colnames(ATAC_mat),"encbias"])
summary_W <- c(DNase_bias[colnames(ATAC_mat),"encbias"])

m <- nrow(Amat)
n <- ncol(Amat)

Amat_sumi <- apply(Amat,2,sum)
Amat_sumj <- apply(Amat,1,sum)
Dmat_sumi <- apply(Dmat,2,sum)
Dmat_sumj <- apply(Dmat,1,sum)

cb_sumj = Amat_sumj + Dmat_sumj

log_this_v <- ATAC_bias[colnames(ATAC_mat),"encbias"]
this_V <- exp((log_this_v - mean(log_this_v)))

log_this_w <- DNase_bias[colnames(ATAC_mat),"encbias"]
this_W <- exp((log_this_w - mean(log_this_w)))

#this_VW <- cbind(this_V,this_W)
#this_V <- exp(ATAC_bias[,"encbias"])#DNase_counts
for(iter_num in 1:maxit){

    ptm <- proc.time()
    ### V,W to U
    this_U <- c()
    for(i in 1:m){
        Ui = cb_sumj[i] / sum((this_V + this_W) * Lmat[i,])
        this_U <- append(this_U,Ui)
    }

    ### U to V,W
    this_V <- c()
    this_W <- c()
    for(j in 1:n){
        Vj = Amat_sumi[j] / sum(this_U * Lmat[,j])
        this_V <- append(this_V,Vj)
        Wj = Dmat_sumi[j] / sum(this_U * Lmat[,j])
        this_W <- append(this_W,Wj)
    }

    summary_U <- cbind(summary_U, this_U)
    summary_V <- cbind(summary_V, this_V)
    summary_W <- cbind(summary_W, this_W)

    #raw_this_V <- this_V
    #log_this_v <- log(raw_this_V)
    #this_V <- exp((log_this_v - mean(log_this_v)))

    #raw_this_W <- this_W
    #log_this_w <- log(raw_this_W)
    #this_W <- exp((log_this_w - mean(log_this_w)))
    print(proc.time() - ptm)

}

rownames(summary_V) <- colnames(ATAC_mat)
colnames(summary_V) <- c("init",paste0("it",seq(maxit)))
rownames(summary_W) <- colnames(ATAC_mat)
colnames(summary_W) <- c("init",paste0("it",seq(maxit)))
rownames(summary_U) <- rownames(ATAC_mat)
colnames(summary_U) <- c("ATAC","DNase",paste0("it",seq(maxit)))

write.table(summary_V,file=paste0(tissue,"_mergePeak_flank",flankN,"_combine_Vstart_Vmat.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(summary_W,file=paste0(tissue,"_mergePeak_flank",flankN,"_combine_Vstart_Wmat.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(summary_U,file=paste0(tissue,"_mergePeak_flank",flankN,"_combine_Vstart_Umat.txt"),row.names=T,col.names=T,quote=F,sep="\t")






#alpha <- rep(seq(M),N)
#beta <- rep(seq(N),each=M)
#
##reg_model <- glm(X~as.factor(alpha)+as.factor(beta),family="poisson")
#reg_model <- glm(X~as.factor(alpha)+as.factor(beta)+L,family="poisson")
##reg_model <- lm(X~as.factor(alpha)+as.factor(beta))
#
#reg_coeff <- reg_model$coefficients
#### 
#Intercept <- reg_coeff["(Intercept)"]
#coeffL <- reg_coeff["L"]
#coeffA <- c(0,reg_coeff[2:M])
#coeffB <- c(0,reg_coeff[(M+1):(M+N-1)])
#
#logV <- c( -1*sum(coeffB)/N )
#for(j in 2:N){
#    logV <- c(logV,coeffB[j]+logV[1])
#}
#
#logU <- c( Intercept + coeffL*L[1] - logV[1])
#for(i in 2:M){
#    logU <- c(logU,coeffA[i]+logU[1])
#}

#ATAC_counts <- apply(ATAC_mat,1,sum)
#
##colnames(X)
#names(logV) <- colnames(ATAC_mat)
#
#write.table(x=logV,file=paste0(tissue,"_mergePeak_flank",flankN,"_",datatype,"_biasV.txt"),row.names=T,col.names=F,sep="\t",quote=F)
#write.table(x=cbind(ATAC_counts, logU),file=paste0(tissue,"_mergePeak_flank",flankN,"_",datatype,"_countU.txt"),row.names=F,col.names=F,sep="\t",quote=F)

#test <- cbind(X,alpha,beta,L)
#U <- matrix(rep(0,M*N*M), nrow=M*N)






