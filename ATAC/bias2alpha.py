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
import copy,time
import numpy
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------


def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
  #      print seq
        nmer_seq[seq] = 0
    del allseq
    return nmer_seq
    
def calculateMN(rawbiasFile):
    Mgenomepos = 0
    Nreadscount = 0
    inf = open(rawbiasFile)
    for line in inf:
        ll = line.split()
        Mgenomepos += int(ll[3])
        Nreadscount += int(ll[2])
    return Nreadscount, Mgenomepos


def generate_alpha(rawbiasF,encbiasF,outfile):
    N,M = calculate(rawbiasF)
    inf = open(encbiasF)
    outf=  open(outfile,'w')
    for line in inf:
        ll = line.split()
        seqtype = ll[0]
        if ll[1] == "NA":
            rawbias = "NA"
            rawalpha = "NA"
        else:
            rawbias = float(ll[1])
            rawaplha = N*1.0/(M*pow(numpy.e,rawbias))
        encbias = float(ll[2])
        encalpha = N*1.0/(M*pow(numpy.e,encbias))
        newll = [seqtype,rawbias,encbias,rawalpha,encalpha]
        outf.write("\t".join(map(str,newll))+"\n")
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
    optparser.add_option("-p","--peak",dest="peak",type="str",
                         help="peak region for calculate background (k-mer occurancy)")
    optparser.add_option("-t","--tag",dest="tag",type="str",
                         help="6-column reads file (bed) for calculate foreground (cuts occyrancy)")              
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    peak = options.peak
    tag = options.tag
    out = options.out
    seq = options.sequence
    flank = options.flank
    if not tag or not peak:
        optparser.print_help()
        sys.exit(1)
    
    seqbias(peak,tag,out,seq,flank)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



