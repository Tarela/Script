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

def GCtagdis(peakfile,taga,tagb,sequence,chromLen,out,delta):
    t = time.time()
    chromDict = readChrom(chromLen)
    pChr,pStart,pEnd,usedChr = readRegion(peakfile,chromDict)
    tagA = readTag(taga,usedChr)
    tagB = readTag(tagb,usedChr)
    genome = twobitreader.TwoBitFile(sequence) 
    #print time.time()-t , "read data"
    t = time.time()
    outf = open(out,'w')
    for i,elem in enumerate(pChr):
        Chr = pChr[i]
        Start = pStart[i]
        End = pEnd[i]
        Seq = genome[Chr]
        if End > len(Seq):
            continue
        A_cuts = tagA[Chr][Start:End][delta:(End-Start-delta)]
        B_cuts = tagB[Chr][Start-1:End-1]   ### here -1 mean transfer -cut to +cut , then the distance is between cuts not tags
        GCcounts = []
        ATcounts = []
        GCalls = []
        ATalls = []
        for i in range(-delta,delta+1):
            if i == 0 : 
                continue
            B_use = B_cuts[(delta+i):(len(B_cuts)-delta+i)]
            if i<0:
                related_loci = [Start+delta+i,End-delta]
            else:
                related_loci = [Start+delta,End-delta+i]
            use_seq = Seq[related_loci[0]:related_loci[1]]
            GCvector = makeGC_numpy(use_seq,abs(i))
            GCcount = (GCvector & A_cuts & B_use).count()
            GCcounts.append(GCcount)
            GCall = GCvector.count()
            ATvector = GCvector
            ATvector.invert()        
            ATcount = (ATvector & A_cuts & B_use).count()
#print GCcount == ATcount
            ATcounts.append(ATcount)
            ATall = ATvector.count()
            GCalls.append(GCall)
            ATalls.append(ATall)
            
        newll = [Chr,Start,End] + GCcounts  + ATcounts + GCalls  + ATalls        
        
        outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
   # print time.time()-t,"process"


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-p","--peakfile",dest="peakfile",type="str",
                         help="")
    optparser.add_option("-a","--taga",dest="taga",type="str",
                         help="")
    optparser.add_option("-b","--tagb",dest="tagb",type="str",
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")
    optparser.add_option("-d","--distance",dest="distance",type="int",default = 100,
                         help="max distance for calculation from +cut to all nearby -cut")

#========minor options=============
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/home/sh430/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-l","--chromLen",dest="chromLen",type="str",default='/mnt/data/static_libraries/chromLen/hg19.len',
                         help="")



    (options,args) = optparser.parse_args()

    peakfile = options.peakfile
    taga = options.taga
    tagb = options.tagb
    delta = options.distance
    sequence = options.sequence
    chromLen = options.chromLen
    out = options.out
    if not peakfile:
        optparser.print_help()
        sys.exit(1)
    
    GCtagdis(peakfile,taga,tagb,sequence,chromLen,out,delta)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

