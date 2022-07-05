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

def addGC(peakfile,sequence,out):

    genome = twobitreader.TwoBitFile(sequence) 
    outf = open(out,'w')
    inf = open(peakfile)
    for line in inf:
        ll = line.split()
        Seq = genome[ll[0]][int(ll[1]):int(ll[2])]
        GC=0
        for i in Seq:
            if i.upper()=="G" or i.upper()=="C":
                GC +=1
        newll = ll + [str(GC)]
        outf.write("\t".join(newll)+"\n")
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
    optparser.add_option("-p","--peakfile",dest="peakfile",type="str",
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="input ll + GC number")
#========minor options=============
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="whole genome sequence in 2bit format")



    (options,args) = optparser.parse_args()

    peakfile = options.peakfile
    sequence = options.sequence
    out = options.out
    if not peakfile:
        optparser.print_help()
        sys.exit(1)
    
    addGC(peakfile,sequence,out)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

