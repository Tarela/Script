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
import twobitreader
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
            r = 'N' #print i 
        revseq += r
    return revseq

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
    
def seqbias(sequence,chrm):

    genome = twobitreader.TwoBitFile(sequence)
    # keep count of the number of occurrences of each n-mer
    
#    wholeSeq = genome[chrm][:]
    for chrm in  genome.keys():
        if not chrm.split('_')[0] ==chrm :
            continue
        print chrm
        wholeSeq = genome[chrm][:]

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
  
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-c","--chrm",dest="chrm",type="str",
                         help="chromsome selected")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    seq = options.sequence
        
    seqbias(seq,options.chrm)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


