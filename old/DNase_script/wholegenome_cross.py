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
import numpy
import bitarray
from math import *


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
        if not tagbit.has_key(ll[0]) and int(ll[2])<= len(tagbit[ll[0]]):
            continue
        if len(ll)>=6 and ll[5] == '-' :
            tagbit[ll[0]][int(ll[2])-1] = 1
        else:
            tagbit[ll[0]][int(ll[1])] = 1
    inf.close()
    return tagbit

def GCtagdis(peakfile,taga,tagb,species,chromLen,out,chrm):
    print 1
    chromDict = readChrom(chromLen)
    print 2
    pChr,pStart,pEnd,usedChr = readRegion(peakfile,chromDict)
    print 3
    tagA = readTag(taga,usedChr)
    print 4
    tagB = readTag(tagb,usedChr)

    cA = tagA[chrm][300:(len(tagA[chrm])-301)]
    cB = tagB[chrm]
    #cA = tagA['chr22'][300:(len(tagA['chr22'])-301)]
    #cB = tagB['chr22']
    #print cA.count()
    #result = []
    
    outf = open(out,'w')
    for i in range(601):
       #print len(cA)
       # print len(cB[i:(len(cB)-601+i)])
     #   result.append((cB[i:(len(cB)-601+i)] & cA).count())
        ct = (cB[i:(len(cB)-601+i)] & cA).count()
        outf.write(str(i-300)+'\t'+ str(ct) +'\n')
    outf.close()
#    print tagA.keys()
#    print len(tagA['chr22'])
#    print tagA['chr22'].count()






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


#========minor options=============
    optparser.add_option("-s","--species",dest="species",type="str",default='hg19',
                         help="")
    optparser.add_option("-l","--chromLen",dest="chromLen",type="str",default='/mnt/data/static_libraries/chromLen/hg19.len',
                         help="")
    optparser.add_option("-c","--chromtain",dest="chromtain",type="str",default='chr22',
                         help="")


    (options,args) = optparser.parse_args()

    peakfile = options.peakfile
    taga = options.taga
    tagb = options.tagb
    species = options.species
    chromLen = options.chromLen
    out = options.out
    if not peakfile:
        optparser.print_help()
        sys.exit(1)
    
    GCtagdis(peakfile,taga,tagb,species,chromLen,out,options.chromtain)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

