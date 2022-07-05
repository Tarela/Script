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
from CistromeAP.taolib.CoreLib.BasicStat.Func import *
from CistromeAP.jianlib.BwReader import BwIO
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
        #print len(v)
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
abline(v=xmax/2+1,lwd=2,lty=2)
legend('topleft',c('plus','minus'),col=c('darkred','darkblue'),lwd=3,lty=c(1,1),bty="n")
legend('topright',c('plusBG','minusBG'),col=c('red','blue'),lwd=3,lty=c(2,2),bty="n")

axis(side=1,at = seq(xmax/2-9,xmax/2+11,10),labels=seq(-10,10,10))
box()
dev.off()
'''%(str(plus)[1:-1],str(minus)[1:-1],str(Pbg)[1:-1],str(Mbg)[1:-1],str(max(plus+minus)),str(min(plus+minus)),len(plus),str(max(Pbg+Mbg)),str(min(Pbg+Mbg)),out+'.pdf')
    Rscript.write(Rstring)
    Rscript.close()
    os.system('Rscript %s'%(out+'.r'))

def make_cut(cutbw_handle,ll,span,flength):
    '''
    flength is fetch length, treat different for +/- strand motif
    return a list
    '''
    total = list(cutbw_handle.summarize(ll[0],max(0,int(ll[1])-flength),int(ll[2])+flength,int(ll[2])-int(ll[1])+ flength *2).sum_data)
    #ts = sum(total)
    if ll[5] != '-':
        out = total[ ( flength - span ) : ( flength + int(ll[2]) - int(ll[1]) + span ) ]
    else:
        out = total[ ( flength - span ) : ( flength + int(ll[2]) - int(ll[1]) + span ) ][::-1]
    return out
    
def readBG(bgmatrix):
    pBG = {}
    nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        ll = line.split()
        name = ll[0]
        pBG[name] = float(ll[1])
        nBG[name] = float(ll[2])
    inf.close()
    return pBG,nBG
    
def getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,fetch_length=100,gen='hg19'):

    
    p=BwIO(pcut)
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    inf = open(inputfile)    
    pp=[]
    pm=[]
    X = c.interval(genome=gen)
    X.chrom,X.start,X.end,X.val = [],[],[],[]
    pBG,nBG = readBG(BGmatrix)
    for line in inf:
        ll = line.split()
        if not chrom_len.has_key(ll[0]):
            continue
        pout = make_cut(pcutbw,ll,pspan,fetch_length)
        nout = make_cut(ncutbw,ll,pspan,fetch_length)
        if ll[5] == "-":
            pout,nout = nout,pout
        if pout == 'NA':
            continue
        #print len(pout),len(nout),ll[:3]
        pp.append(pout)
        pm.append(nout)
        X.chrom.append(ll[0])
        X.start.append(int(ll[1])-pspan -3   + 1)
        X.end.append(int(ll[2]) + pspan +3   + 1)
        X.val.append(ll[5])
#total[ ( flength - span ) : ( flength + int(ll[2]) - int(ll[1]) + span ) ]

    meanp = apply_mean(pp)
    meanm = apply_mean(pm) 

    X.getSequence()
    
    pbglist = []
    nbglist = []
    for i,elem in  enumerate(X.seq):
        seq = X.seq[i]
        strand = X.val[i]
        if 'N' in seq.upper():
            continue
        pseq = seq[:-1]
        nseq = seq[1:]
        #if 'N' in pseq  or 'N' in nseq:
        #    continue
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - 6):
            p.append(pBG[pseq[k:k+6].upper()])
            n.append(nBG[nseq[k:k+6].upper()])
        if strand != '-':
            pbglist.append(p)
            nbglist.append(n)
        else:
            pbglist.append(n[::-1])
            nbglist.append(p[::-1])
    #print nbglist
    meanpbglist = apply_mean(pbglist)
    meanmbglist = apply_mean(nbglist)        

    plot_template(meanp,meanm,meanpbglist,meanmbglist,outputfile)



# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")                         
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="")
    optparser.add_option("-p","--positive_cut",dest="pcut",type="str",
                         help="")
    optparser.add_option("-n","--negative_cut",dest="ncut",type="str",
                         help="")

    optparser.add_option("--plot_span",dest="pspan",type="int",default = 10,
                         help="default = 10 , means total length = 10+bedlen+10")

    optparser.add_option("--genome",dest="genome",type="str",default = 'hg19',
                         help="choose from hg19 and mm9")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    bgmatrix = options.bgmatrix
    pcut = options.pcut
    ncut = options.ncut
    genome = options.genome
#    method = options.method
    pspan = options.pspan
    
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,bgmatrix,pcut,ncut,pspan,fetch_length=100,gen=genome)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
