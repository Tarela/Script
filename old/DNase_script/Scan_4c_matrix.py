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
def makeTCFOS(cutbw_plus,cutbw_minus,ll,tspan,ml):
    
    pcut = list(cutbw_plus.summarize(ll[0], ( int(ll[1]) + int(ll[2]) )/2 -tspan ,( int(ll[1]) + int(ll[2]) )/2 +tspan ,2*tspan).sum_data)
    ncut = list(cutbw_minus.summarize(ll[0], ( int(ll[1]) + int(ll[2]) )/2 -tspan ,( int(ll[1]) + int(ll[2]) )/2 +tspan ,2*tspan).sum_data)
    TC =  sum(pcut)+sum(ncut)

    C = sum(pcut[(tspan-ml/2) : (tspan-ml/2+ml)]) + sum(ncut[(tspan-ml/2) : (tspan-ml/2+ml)])
    L = sum(pcut[(tspan-ml/2-ml):(tspan-ml/2)]) + sum(ncut[(tspan-ml/2-ml):(tspan-ml/2)])
    R = sum(pcut[(tspan-ml/2+ml):(tspan-ml/2+2*ml)]) + sum(ncut[(tspan-ml/2+ml):(tspan-ml/2+2*ml)])
 
    FOS = ( (C+1)/(R+1) + (C+1)/(L+1) )
    return TC,FOS

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
    
def getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,tspan,gen,fetch_length=100):

    
 #   p=BwIO(pcut)
 #   chrom_len = {}
 #   for i in p.chromosomeTree['nodes']:
 #       chrom_len[i['key']] = i['chromSize']
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
        seq = genome[chrom][(start-pspan-3):(end + pspan+3)]
        pout = make_cut(pcutbw,ll,pspan,fetch_length)
        nout = make_cut(ncutbw,ll,pspan,fetch_length)
    
        if strand == "-":
            pout,nout = nout,pout
            
        if pout == 'NA':
            continue        

        if 'N' in seq.upper():
            continue
        #print 1
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - 6):
            p.append(pBG[pseq[k:k+6].upper()])
            n.append(nBG[nseq[k:k+6].upper()])
        if strand != '-':
            pbglist = p
            nbglist = n
        else:
            pbglist = n[::-1]
            nbglist = p[::-1]
        TC,FOS = makeTCFOS(pcutbw,ncutbw,ll,tspan,ml)
        newll = ll  + [TC,FOS] + pout + nout  + pbglist + nbglist
        outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
    inf.close()
 
    


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

    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 25,
                         help="default = 25 , means total length = 50")
    optparser.add_option("-t","--tagcount_span",dest="tspan",type="int",default = 100,
                         help="default = 100 , means total length = 200")

    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="2bit format")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    bgmatrix = options.bgmatrix
    pcut = options.pcut
    ncut = options.ncut
    genome = options.genome
    tspan = options.tspan
#    method = options.method
    pspan = options.pspan
  
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,bgmatrix,pcut,ncut,pspan,tspan,genome,fetch_length=100)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
