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
import math,time
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader

import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def rev(seq):
    revseq = ""
    for i in seq[::-1]:
        if i == 'A':
            r = 'T'
        elif i == 'T':
            r = 'A'
        elif i == 'C':
            r = 'G'
        elif i == 'G':
            r = 'C'
        else:
            print i
        revseq += r
    return revseq

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
    
def readChromLen_wig(chromLen,bwtrack,usechrom):
    bwhandle = BigWigFile(open(bwtrack, 'rb'))
    chromLen_dict = {}
    inf = open(chromLen)
    for line in inf:
        ll = line.split()
        if ll[0] == usechrom or usechrom == "all":
            track = list(bwhandle.summarize(ll[0],0,int(ll[1]),int(ll[1])).sum_data)
            chromLen_dict[ll[0]] = track
    inf.close()
    return chromLen_dict
    
    #### now here

def sitepro_scan(out,w_plus,w_minus,chromLen,bgmatrix,span,gen,lflank,rflank,usechrom):
    bwdict_plus = readChromLen_wig(chromLen,w_plus,usechrom)
    bwdict_minus = readChromLen_wig(chromLen,w_minus,usechrom)
    for i in bwdict_plus.keys():
        print i, len(bwdict_plus[i])
    #sys.exit(1)
    nmer = lflank + rflank
    genome = twobitreader.TwoBitFile(gen)
    pBG,nBG = readBG(bgmatrix)
    outf = open(out,'w')
    for chrm in bwdict_plus.keys():
        allcuts = bwdict_plus[chrm]
        for bp in range(len(allcuts)):
            cut = allcuts[bp]
            seq = genome[chrm][(bp-lflank-span):(bp+rflank+span)]
            pseq = seq[:-1]
            window_bias = []
            for i in range(len(pseq) + 1 -lflank-rflank):
                window_bias.append( pBG[pseq[i:(i+6)].upper()] )
            window_count = sum(bwdict_plus[chrm][(bp-span):(bp+span)])
            if len(window_bias) == 2*span:
                predict_cut = window_count*1.0*window_bias[span]/sum(window_bias)
                newll = [cut,window_count,window_bias[span],predict_cut]
                outf.write("\t".join(map(str,newll))+"\n")
            else:
                print 'plus',(bp-lflank-span),(bp+rflank+span),len(genome[chrm][:]),bp,len(window_bias),len(pseq)
    for chrm in bwdict_minus.keys():
        allcuts = bwdict_minus[chrm]
        for bp in range(len(allcuts)):
            cut = allcuts[bp]
            seq = genome[chrm][(bp-lflank-span):(bp+rflank+span)]
            nseq = seq[1:]
            window_bias = []
            for i in range(len(nseq) + 1 -lflank-rflank):
                window_bias.append( nBG[nseq[i:(i+6)].upper()] )
            window_count = sum(bwdict_plus[chrm][(bp-span):(bp+span)])
            if len(window_bias) == 2*span:
                predict_cut = window_count*1.0*window_bias[span]/sum(window_bias)
                newll = [cut,window_count,window_bias[span],predict_cut]
                outf.write("\t".join(map(str,newll))+"\n")
            else:
                print 'minus',(bp-lflank-span),(bp+rflank+span),len(genome[chrm][:]),bp,len(window_bias),len(nseq)

    ### now here
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

    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/DNase_paper/Data/Yeast_DNase/yeast_naked_plus.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/DNase_paper/Data/Yeast_DNase/yeast_naked_minus.bw",
                         help="")
    optparser.add_option("-l","--chromLen",dest="chromLen",type="str",default = "/home/sh430/Data/Genome/ChromLen/ChromInfo_Scerevisiae.txt",
                         help="chromsome length file")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",default = "/home/sh430/Project/DNase_paper/Data/BGmatrix_new/yeast_naked_wholegenome_newbg.txt",
                         help="sequence bias matrix")
#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 25,
                         help="region for get total signal in single bp, default = 25 means +-25bp(total 50bp) signal as total for each bp")
    optparser.add_option("--genome",dest="genome",type="str",default = "/home/sh430/Data/Genome/Scerevisiae.2bit",
                         help="2bit format")
    optparser.add_option("--left",dest="leftflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--right",dest="rightflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--usechrom",dest="usechrom",type="str",default = "all",
                         help="")


    (options,args) = optparser.parse_args()

    out=  options.out
    w_plus = options.w_plus
    w_minus = options.w_minus
    chromLen = options.chromLen
    bgmatrix = options.bgmatrix
    gen = options.genome
    lflank = options.leftflank
    rflank = options.rightflank
    usechrom = options.usechrom
    
    sitepro_scan(out,w_plus,w_minus,chromLen,bgmatrix,options.Cspan,gen,lflank,rflank,usechrom)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

