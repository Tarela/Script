
tissue <- "Ovary"


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
for(datatype in c("ATAC","DNase")){
  pdf(file=paste0(tissue,"_",datatype,"Ceffect.pdf"),width=15,height=10)
  par(mfrow=c(3,5),mar=c(4,4,2,2))

  for(flankN in 1:3){
    for(codetype in c("iterPos","samplingPos","ZIP","simplexZIPgroup","GDsimplex")){
      outU <- read.table(paste0(paste(tissue,flankN,datatype,"outU",codetype,sep="_"),".txt"), row.names=1)
      ATAC_ratio <- as.numeric(outU[,2]/(outU[,1]+1))
      names(ATAC_ratio)<-rownames(outU)
      boxplot(ATAC_ratio,ATAC_ratio[TSSpeak],ATAC_ratio[TTSpeak],ATAC_ratio[otherpeak],outline=F,names=c("all","TSS","TTS","other"),ylab="U/reads",las=2,main=codetype)
      ATTS <- -log10(t.test(ATAC_ratio[TTSpeak],ATAC_ratio,alternative="less")$p.val)
      ATSS <- -log10(t.test(ATAC_ratio[TSSpeak],ATAC_ratio,alternative="less")$p.val)
      Aother <- -log10(t.test(ATAC_ratio[otherpeak],ATAC_ratio,alternative="less")$p.val)
      legend("topleft",legend=paste0(c("TSS","TTS","other"),":",c(ATSS,ATTS,Aother)),bty="n")
    }
  }
  dev.off()
}
