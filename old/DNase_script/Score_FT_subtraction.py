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
import twobitreader
#from CistromeAP.taolib.CoreLib.BasicStat.Func import *
#from CistromeAP.jianlib.BwReader import BwIO
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
    
    
def getsignal(inputfile,outputfile,cut,pcut,ncut,pspan,BGmatrix,gen,left,right):

    genome = twobitreader.TwoBitFile(gen)
#    p=BwIO(pcut)
#    chrom_len = {}
#    for i in p.chromosomeTree['nodes']:
#        chrom_len[i['key']] = i['chromSize']
    cutbw = BigWigFile(open(cut, 'rb'))
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    inf.seek(0)
    outf = open(outputfile,'w')
    pBG,nBG = readBG(BGmatrix)
    for line in inf:
        ll = line.split()

        ## make tag count
        cut = list(cutbw.summarize(ll[0],int(ll[1]) + ml/2 -pspan ,int(ll[1]) + ml/2 +pspan ,2*pspan).sum_data)
        TC = sum(cut)
        ## make L, C , R
        pcut = list(pcutbw.summarize(ll[0],int(ll[1]) - ml ,int(ll[2]) + ml ,3*ml).sum_data)
        ncut = list(ncutbw.summarize(ll[0],int(ll[1]) - ml ,int(ll[2]) + ml ,3*ml).sum_data)
        L = pcut[:ml] + ncut[:ml]
        C = pcut[ml: (ml*2)] + ncut[ml: (ml*2)]
        R = pcut[(ml*2):(ml*3)] + ncut[(ml*2):(ml*3)]
        ## make bias for L,C,R
        seq = genome[ll[0]][(int(ll[1]) - ml - left):(int(ll[2]) + ml + right)]
        if 'N' in seq.upper():
            continue
        #print 1
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - left-right):
            p.append(pBG[pseq[k:k+left+right].upper()])
            n.append(nBG[nseq[k:k+left+right].upper()])
        Lbg = p[:ml] + n[:ml]
        Cbg = p[ml: (ml*2)] + n[ml: (ml*2)]
        Rbg = p[(ml*2):(ml*3)] + n[(ml*2):(ml*3)]
        ## correction , cleavage/bias
        Lcr = []
        Ccr = []
        Rcr = []
        for i in range(len(L)):
            Lcr.append(float(L[i])/float(Lbg[i]))
            Ccr.append(float(C[i])/float(Cbg[i]))
            Rcr.append(float(R[i])/float(Rbg[i]))

        FOS = -1*( (sum(C)+1)/(sum(R)+1) + (sum(C)+1)/(sum(L)+1) ) 
        FOScr = -1*( (sum(Ccr)+1)/(sum(Rcr)+1) + (sum(Ccr)+1)/(sum(Lcr)+1) )
        
        newll = ll + [TC,FOS,FOScr]
        outf.write("\t".join(map(str,newll))+"\n")

    outf.close()
    


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

    optparser.add_option("-g","--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="")

    optparser.add_option("-w","--cutbw",dest="cut",type="str",
                         help="")
    optparser.add_option("-p","--pcutbw",dest="pcut",type="str",
                         help="")
    optparser.add_option("-n","--ncutbw",dest="ncut",type="str",
                         help="")                         
    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 100,
                         help="window size for tagcount and dymDHS, default = 100 , means total length = 200")
    optparser.add_option("--left",dest="left",type="int",default = 3,
                         help="")
    optparser.add_option("--right",dest="right",type="int",default =3,
                         help="")                         

 
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    cutbw = options.cut
    pcutbw = options.pcut
    ncutbw = options.ncut
    BGmatrix = options.bgmatrix
    pspan = options.pspan
    gen = options.genome
    left = options.left
    right = options.right
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,cutbw,pcutbw,ncutbw,pspan,BGmatrix,gen,left,right)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
