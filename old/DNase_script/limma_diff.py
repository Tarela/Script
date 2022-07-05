import os
def makelimma(inputfile,treatnum,ctrlnum,outfile,limmaname):
    outf = open(outfile,'w')
    outf.write("""
library(limma)
data<-read.table("%s",sep="\t",header=T)
exp<-log(data[,3:ncol(data)] + 0.01)
gm.design<-cbind(ctrl=c(%s),treat=c(%s))
gm.fit<-lmFit(exp,gm.design)
gm.matrix<-makeContrasts(CvsL=treat-ctrl,levels=gm.design)
gm.fit<-contrasts.fit(gm.fit,gm.matrix)
gm.fit <- eBayes(gm.fit)
t<-topTable(gm.fit,adjust="BH",number=40000,p.value=1)
### using limma to make differential expression gene
out<-cbind(data[,1:2],t[,2:ncol(t)])
write.table(out,file="%s",quote=F,sep="\t",row.names=F)
    """%(inputfile,str([0]*treatnum)[1:-1]+','+str([1]*ctrlnum)[1:-1],str([1]*treatnum)[1:-1]+','+str([0]*ctrlnum)[1:-1],limmaname))
    outf.close()
    os.system('Rscript %s'%(outfile))
    return 1
makelimma('Expression.txt',2,2,'limmatext.r','limmatext.txt')


