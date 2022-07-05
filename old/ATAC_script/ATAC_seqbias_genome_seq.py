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
    
def seqbias(out,sequence,rflank,lflank):

    genome = twobitreader.TwoBitFile(sequence)
    # keep count of the number of occurrences of each n-mer
    
    seq_nmer_dict   = make_nmer_dict(lflank+rflank)
    for chrom in genome.keys():
        if not chrom.split('_')[0] == chrom :
            continue
        wholeSeq = genome[chrom][:]
        RVwholeSeq = rev(wholeSeq)
        for i in range(len(wholeSeq)-lflank-rflank +1 ):
            seq6mer = wholeSeq[i:(i+lflank+rflank)]
            RVseq6mer = RVwholeSeq[i:(i+lflank+rflank)]
            if 'a' in seq6mer or 't' in seq6mer or 'c' in seq6mer or 'g' in seq6mer or 'n' in seq6mer or 'N' in seq6mer : 
                pass
            else:
                seq_nmer_dict[seq6mer] += 1
            if 'a' in RVseq6mer or 't' in RVseq6mer or 'c' in RVseq6mer or 'g' in RVseq6mer or 'n' in RVseq6mer or 'N' in RVseq6mer : 
                pass
            else:
                seq_nmer_dict[RVseq6mer] += 1


    outf = open(out,'w')
    for seqtype in sorted(seq_nmer_dict.keys()):
        outf.write("\t".join(map(str,[seqtype,seq_nmer_dict[seqtype]]))+"\n")
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
         
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("--left",dest="left",type="int",default = 3,
                         help="region length(left) for calculate seqbias , default =3  nmer = left+right")
    optparser.add_option("--right",dest="right",type="int",default = 3,
                         help="region length(right) for calculate seqbias , default =3  nmer = left+right")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    out = options.out
    seq = options.sequence
    lflank = options.left
    rflank = options.right
    if not out:
        optparser.print_help()
        sys.exit(1)
    
    seqbias(out,seq,rflank,lflank)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


