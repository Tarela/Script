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
import twobitreader
import numpy
import random
#from CistromeAP.taolib.CoreLib.BasicStat.Func import *
#from CistromeAP.jianlib.BwReader import BwIO
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
#import scipy.stats.distributions

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

def getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,tspan,gen,left,right):
    pBG,nBG = readBG(BGmatrix)

    genome = twobitreader.TwoBitFile(gen)
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))

    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    pspan = pspan - ml/2
    inf.seek(0)
    pBG,nBG = readBG(BGmatrix)
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()

        chrom = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        strand = ll[5]
        
        if int(ll[1])-pspan-tspan < 0 :
            continue
        pout = list(pcutbw.summarize(ll[0],max(0,int(ll[1])-pspan-tspan),int(ll[2])+ pspan + tspan,int(ll[2])-int(ll[1])+ pspan *2 + tspan*2).sum_data)
        nout = list(ncutbw.summarize(ll[0],max(0,int(ll[1])-pspan-tspan),int(ll[2])+ pspan + tspan,int(ll[2])-int(ll[1])+ pspan *2 + tspan*2).sum_data)
        if strand != '-':
            pass
        else:
            pout,nout = nout[::-1],pout[::-1]
        seq = genome[chrom][(start-pspan-tspan-left):(end + pspan + tspan+right)]
        ### get seqbias
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
        if strand != '-':
            pbglist = p
            nbglist = n
        else:
            pbglist = n[::-1]
            nbglist = p[::-1]
        
        if len(pbglist) != len(pout):
            print 'len pbglist != len pout , should have bugs'
            sys.exit(1)
            
        p_cleavage = []
        n_cleavage = []
        p_bias =   []
        n_bias = []
        p_predict = []
        n_predict = []
        p_uniform = []
        n_uniform = []
        
        for i in range(len(pout)-2*tspan):
            p_ob = pout[i+tspan]
            n_ob = nout[i+tspan]
            ptotal = sum(pout[i:(i+tspan*2)])
            ntotal = sum(nout[i:(i+tspan*2)])
            pbias = pbglist[i+tspan]
            nbias = nbglist[i+tspan]
            pbias_total = sum(pbglist[i:(i+tspan*2)])
            nbias_total = sum(nbglist[i:(i+tspan*2)])
            pbias_por = pbias*1.0/pbias_total
            nbias_por = nbias*1.0/nbias_total
            ppredict = pbias_por * ptotal
            npredict = nbias_por * ntotal
            puniform = ptotal*1.0/(2*tspan)
            nuniform = ntotal*1.0/(2*tspan)

            p_cleavage.append(p_ob)
            n_cleavage.append(n_ob)
            p_bias.append(pbias)
            n_bias.append(nbias)
            p_predict.append(ppredict)
            n_predict.append(npredict)
            p_uniform.append(puniform)
            n_uniform.append(nuniform)

        newll = ll + p_cleavage +  n_cleavage + p_bias +  n_bias+  p_uniform + n_uniform + p_predict + n_predict 
        outf.write("\t".join(map(str,newll)) + "\n")
    inf.close()
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
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="")
    optparser.add_option("-p","--positive_cut",dest="pcut",type="str",
                         help="")
    optparser.add_option("-n","--negative_cut",dest="ncut",type="str",
                         help="")
    optparser.add_option("--left",dest="left",type="int",default=3,
                         help="")
    optparser.add_option("--right",dest="right",type="int",default=3,
                         help="")


    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 25,
                         help="default = 25 , means total plot span length = 50")
    optparser.add_option("-t","--tag_span",dest="tspan",type="int",default = 25,
                         help="default = 25 , means total tag span length = 50")

    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="2bit format")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    BGmatrix = options.bgmatrix
    pcut = options.pcut
    ncut = options.ncut
    genome = options.genome
    pspan = options.pspan
    tspan = options.tspan

    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
#    getsignal(inputfile,outputfile,bgmatrix,pcut,ncut,pspan,genome,options.left,options.right,fetch_length=100)
    #getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,genome,options.left,options.right,matrix_p,matrix_n,options.prelim,options.cutlim,fetch_length=100)
    getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,tspan,genome,options.left,options.right)
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
