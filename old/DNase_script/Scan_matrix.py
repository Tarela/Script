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

def make_cut(cutbw_handle,ll,span,flength):
    '''
    flength is fetch length, treat different for +/- strand motif
    return a list
    '''
    total = list(cutbw_handle.summarize(ll[0],max(0,int(ll[1])-flength),int(ll[2])+flength,int(ll[2])-int(ll[1])+ flength *2).sum_data)
    #ts = sum(total)
    if ll[3] != '-':
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

    
 #   p=BwIO(pcut)
 #   chrom_len = {}
 #   for i in p.chromosomeTree['nodes']:
 #       chrom_len[i['key']] = i['chromSize']
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    pspan = pspan - ml/2
    inf.seek(0)
    X = c.interval(genome=gen)
    X.chrom,X.start,X.end,X.val = [],[],[],[]
    pBG,nBG = readBG(BGmatrix)
    for line in inf:
        ll = line.split()
 #       if not chrom_len.has_key(ll[0]):
 #           continue
        
        X.chrom.append(ll[0])
        X.start.append(int(ll[1])-pspan -3   + 1)
        X.end.append(int(ll[2]) + pspan +3   + 1)
        X.val.append(ll[5])
    inf.close()
    X.getSequence()
    
    outf = open(outputfile,'w')
    for i,elem in  enumerate(X.seq):
        
        pchrm = X.chrom[i]
        pstart = X.start[i] -1 +3 + pspan
        pend = X.end[i] -1 -3 - pspan
        seq = X.seq[i]
        strand = X.val[i]
        pll = [pchrm,pstart,pend,strand]
        pout = make_cut(pcutbw,pll,pspan,fetch_length)
        nout = make_cut(ncutbw,pll,pspan,fetch_length)
        if pll[3] == "-":
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
    #print nbglist
        newll = [pchrm,pstart,pend,strand] + pout + nout + pbglist + nbglist
        #print len(pout),len(nout),len(pbglist),len(nbglist),len(newll)
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
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="")
    optparser.add_option("-p","--positive_cut",dest="pcut",type="str",
                         help="")
    optparser.add_option("-n","--negative_cut",dest="ncut",type="str",
                         help="")

    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 25,
                         help="default = 25 , means total length = 50")

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
