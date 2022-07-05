#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re,time
from optparse import OptionParser
import logging
import string
import numpy
import bitarray
from math import *
import twobitreader
import cProfile
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def readChrom(chromLen):
    chromDict = {}
    inf = open(chromLen)
    for line in inf:
        ll = line.split()
        if  len(ll) == 2:
            chromDict[ll[0]] = int(ll[1])
    inf.close()
    return chromDict

def readRegion(peakfile,chromDict): ## keep the peak merged
    usedChr = {}
    usedChrout = {}
    peakChr = []
    peakStart = []
    peakEnd = []
    inf = open(peakfile)
    for line in inf:
        ll = line.split()
        peakChr.append(ll[0])
        usedChr[ll[0]] = 1
        peakStart.append(int(ll[1]))
        peakEnd.append(int(ll[2]))
    inf.close()
    for chrm in usedChr.keys():
        if chromDict.has_key(chrm):
            usedChrout[chrm] = chromDict[chrm]
    return peakChr,peakStart,peakEnd,usedChrout
    
def readTag(tagfile,usedChr):
    inf = open(tagfile)
    tagbit = {}
    for chrm in usedChr.keys():
        chrmbit = bitarray.bitarray(usedChr[chrm])
        chrmbit.setall(False)
        tagbit[chrm] = chrmbit
    for line in inf:
        ll = line.split()
        if len(ll) == 0:
            continue
        if not tagbit.has_key(ll[0]): 
            continue
        if not int(ll[2])<= len(tagbit[ll[0]]):
            continue
        if len(ll)>=6 and ll[5] == '-' :
            tagbit[ll[0]][int(ll[2])-1] = 1
        else:
            tagbit[ll[0]][int(ll[1])] = 1
    inf.close()
    return tagbit

def makeGC(use_seq,length):
    cutoff = length*1.0/2
    GCvector = bitarray.bitarray()
    GCs= []
   # print len(use_seq)
   # print length
    #print len(range(len(use_seq)-length))
    #print len(use_seq)
    for bp in range(len(use_seq)-length):
        seq = use_seq[bp:(bp+length)]
        GCthis = 0
        for i in seq:
            if i.upper() in ['G','C']:
                GCthis += 1
        GCs.append(GCthis)
        if GCthis >= cutoff :
            GCvector.append(1)
        else:
            GCvector.append(0)
    return GCvector
    
def makeGC_numpy(use_seq,length):
    cutoff = length*1.0/2
    GCvector = bitarray.bitarray()
    seq = numpy.array(list(use_seq.upper()))
    b = (seq == 'C')
    d = (seq == 'G')
    f = 1*(b|d)
    cs = numpy.hstack(([0],f.cumsum()))
    GCv = cs[(length):(len(cs)-1)] - cs[:(len(cs)-length-1)]
    GCvector.extend(GCv>=cutoff)
    return GCvector

def GCtagdis(peakfile,out,lowlim,highlim):

    inf = open(peakfile)
    score = {}
    for i in range(lowlim,highlim):
        score[i]=[0]*(i+1)
    for line in inf:
        ll = line.split()
        fraglen = int(ll[3])
        GC = int(ll[4])
        try :
            score[fraglen][GC] += 1
        except:
            pass
    outf1 = open(out,'w')
    newll = ['fraglen'] + range(highlim)
    outf1.write("\t".join(map(str,newll))+"\n")
    for i in range(lowlim,highlim):
        newll = [i]+score[i]
        if len(newll) < 301:
            newll += [0] * (301-len(newll))
        outf1.write("\t".join(map(str,newll))+"\n")

    outf1.close()
 

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-p","--PEtag",dest="PEtag",type="str",
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="raw table ,x=fraglen, y=bp number that is GC, z=count")
#    optparser.add_option("-t","--out2",dest="out2",type="str",
#                         help="percentage table , x=fraglen, y=% of GC in such frag, z=count")
#========minor options=============
    optparser.add_option("--lowlim",dest="lowlim",type="int",default=10,
                         help="")
    optparser.add_option("--highlim",dest="highlim",type="int",default=300,
                         help="")



    (options,args) = optparser.parse_args()

    peakfile = options.PEtag
    out = options.out
 #   out2 = options.out2
    if not peakfile:
        optparser.print_help()
        sys.exit(1)
    
    GCtagdis(peakfile,out,options.lowlim,options.highlim)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


