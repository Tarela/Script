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
from math import *
import twobitreader
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def addGC(inputfile,sequence,out):

    genome = twobitreader.TwoBitFile(sequence) 
    outf = open(out,'w')
    inf = open(inputfile)
    for line in inf:
        ll = line.split()
#        center = (int(ll[1]) + int(ll[2]))/2
        start = int(ll[1])#max(0,center-extsize)
        end = int(ll[2])#center+extsize
        Seq = genome[ll[0]][start:end]
        GC=0
        for i in Seq:
            if i.upper()=="G" or i.upper()=="C":
                GC +=1
        newll = ll + [str(GC*1.0/len(Seq))]
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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="input ll + GC number")
#========minor options=============
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/scratch/sh8tv/Data/Genome/hg38/hg38.2bit',
                         help="whole genome sequence in 2bit format")
    ##optparser.add_option("-e","--extsize",dest="extsize",type="int",default=500,
     #                    help="extsize from the center, default is 500")



    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    sequence = options.sequence
    out = options.out
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    addGC(inputfile,sequence,out)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

