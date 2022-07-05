'''
Created on XXXX-XX-XX

@author: Tarela
'''
#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string
import math
from cistrome import regions as c
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def apply_mean(vector):
    mv = [0]*len(vector[0])
    for v in vector:
        for i in range(len(mv)):
            mv[i] += float(v[i])
    for i in range(len(mv)):
        mv[i] = mv[i]*1.0/len(vector)

    return mv
def plot_template(plus,minus,Pbg,Mbg,out):
    Rscript = open(out+'.r','w')
    Rstring = '''
plus<-c(%s)
minus<-c(%s)
pbg<-c(%s)
mbg<-c(%s)
ymax<-%s
ymin<-%s
xmax<-%s
ybgmax<-%s
ybgmin<-%s
pdf(file = "%s")
plot(plus,type="l",col="darkred",xlim=c(0,xmax),ylim=c(ymin,ymax),lwd=2,xlab="Relative distance (bp)",ylab="Average profile",main="",axes=F)
par(new=T)
plot(minus,type="l",col="darkblue",xlim=c(0,xmax),ylim=c(ymin,ymax),lwd=2,xlab="",ylab="",main="",axes=F)
axis(side=2)
par(new=T)
plot(pbg,type="l",col="red",lty=2,xlim=c(0,xmax),ylim=c(ybgmin,ybgmax),lwd=2,xlab="",ylab="",main="",axes=F)
par(new=T)
plot(mbg,type="l",col="blue",lty=2,xlim=c(0,xmax),ylim=c(ybgmin,ybgmax),lwd=2,xlab="",ylab="",main="",axes=F)
axis(side=4)
abline(v=xmax/2,lwd=2,lty=2)
legend('topleft',c('plus','minus','plusBG','minusBG'),col=c('darkred','darkblue','red','blue'),lwd=3,lty=c(1,1,2,2),bty="n")
axis(side=1,at = seq(0,xmax,10),labels=seq(-xmax/2,xmax/2,10))
box()
dev.off()
'''%(str(plus)[1:-1],str(minus)[1:-1],str(Pbg)[1:-1],str(Mbg)[1:-1],str(max(plus+minus)),str(min(plus+minus)),len(plus),str(max(Pbg+Mbg)),str(min(Pbg+Mbg)),out+'.pdf')
    Rscript.write(Rstring)
    Rscript.close()
    os.system('Rscript %s'%(out+'.r'))


def make_template(data,flank,pflank,topmotif,out,pbw,mbw,bgmatrix,gen):
    w_plus_H=BigWigFile(open(pbw, 'rb'))
    w_minus_H=BigWigFile(open(mbw, 'rb'))
    i =0
    templatelist = []
    pp=[]
    pm=[]
    inf = open(data)
    l1st = inf.readline().split()
    ml = int(l1st[2])-int(l1st[1])
    inf.seek(0)
    for line in inf:
        #if i >= topmotif:
         #   break
        ll = line.split()
        templatelist.append(ll)

    inf.close()
    templatelist.sort(key = lambda x:float(x[4]),reverse=True)

### for cut sitepro
    for ll in templatelist:
        p_sum = list(w_plus_H.summarize(ll[0],int(ll[1])-flank,int(ll[1])+flank,2*flank).sum_data)
        m_sum = list(w_minus_H.summarize(ll[0],int(ll[1])-flank,int(ll[1])+flank,2*flank).sum_data)
        if ll[5] == "+":
            pp.append(p_sum[(flank + 1 + ml/2 -pflank):(flank +1 + ml/2  + pflank)])
            pm.append(m_sum[(flank +1 + ml/2 - pflank):(flank +1 + ml/2  + pflank)])
        if ll[5] == '-' :
            pm.append(p_sum[::-1][(flank +1 + ml/2 -1 - ml - pflank) : (flank +1 + ml/2 -1 -ml +pflank)])
            pp.append(m_sum[::-1][(flank +1 + ml/2 -1 - ml - pflank) : (flank +1 + ml/2 -1 -ml +pflank)])
    print pp
    print pm
    meanp = apply_mean(pp)
    meanm = apply_mean(pm)
    allsum = sum(meanp)+sum(meanm)
    P=[]
    M=[]
    for i in range(len(meanp)):
        P.append(meanp[i])#/allsum)
        M.append(meanm[i])#/allsum)    

### for seqbias bg
    pBG = {}
    nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        ll = line.split()
        name = ll[0]
        pBG[name] = float(ll[1])
        nBG[name] = float(ll[2])
    inf.close()
    X = c.interval(genome=gen)
    X.chrom,X.start,X.end,X.val = [],[],[],[]
    for ll in templatelist:
        X.chrom.append(ll[0])
        X.start.append(int(ll[1])+1-flank)
        X.end.append(int(ll[1])+1+flank)
        X.val.append(ll[5])
    X.getSequence()
    
    pbglist = []
    nbglist = []
    for i,elem in  enumerate(X.seq):
        seq = X.seq[i]
        strand = X.val[i]
        if strand != '+' or 'N' in seq or 'n' in seq:
            continue
        pseq = seq[(flank + 1 + ml/2 -pflank -3):(flank +1 + ml/2  + pflank +2)]
        nseq = seq[(flank + 1 + ml/2 -pflank -2):(flank +1 + ml/2  + pflank +3)]
        #if 'N' in pseq  or 'N' in nseq:
        #    continue
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - 6):
            p.append(pBG[pseq[k:k+6].upper()])
            n.append(nBG[nseq[k:k+6].upper()])
        pbglist.append(p)
        nbglist.append(n)
    
    #print pbglist
    #print nbglist
    meanpbglist = apply_mean(pbglist)
    meanmbglist = apply_mean(nbglist)
    allsum = sum(meanpbglist)+sum(meanmbglist)
    Plusbg=[]
    Minusbg=[]
    for i in range(len(meanpbglist)):
        Plusbg.append(meanpbglist[i])#/allsum)
        Minusbg.append(meanmbglist[i])#/allsum)    
    
        
    plot_template(P,M,Plusbg,Minusbg,out)
   # plot_template(P,M,P,M,out)
   
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputdata",dest="data",type="str",
                         help="")
    optparser.add_option("-o","--outputdata",dest="out",type="str",
                         help="")
    optparser.add_option("-p","--pbw",dest="pbw",type="str",
                         help="")
    optparser.add_option("-m","--mbw",dest="mbw",type="str",
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bg",type="str",
                         help="")
#========minor options=============
    optparser.add_option("--flank",dest="flank",type="int",default = 100,
                         help="+- N for flanking region , default +- 100")
    optparser.add_option("--pflank",dest="pflank",type="int",default = 20,
                         help="+- N for plot flanking region default +- 20")
    optparser.add_option("--genome",dest="genome",type="str",default = 'hg19',
                         help="refgeonme , default hg19 , choice mm9")


    optparser.add_option("--topmotif",dest="topmotif",type="int",default = 5000,
                         help="motif number for making tamplate, default is 5000")


    (options,args) = optparser.parse_args()

    peak = options.data
    out = options.out
    flank = options.flank
    topmotif = options.topmotif
    pflank = options.pflank
    pbw = options.pbw
    mbw = options.mbw
    bg = options.bg
    genome = options.genome
    

    if not peak :
        optparser.print_help()
        sys.exit(1)
    
    make_template(peak,flank,pflank,topmotif,out,pbw,mbw,bg,genome)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


