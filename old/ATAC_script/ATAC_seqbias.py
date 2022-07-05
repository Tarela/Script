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
def readseqbg(seqbg):
    genomeBG = {}
    cutBG = {}
    inf = open(seqbg)
    for line in inf:
        ll = line.split()
        genomeBG[ll[0]] = float(ll[1])
        cutBG[ll[0]] = 0
    inf.close()
    return genomeBG,cutBG
    
def seqbias(tag,out,sequence,lflank,rflank,seqbg):
    genome = twobitreader.TwoBitFile(sequence) 
    genomeBG,cutBG = readseqbg(seqbg)
    
    inf = open(tag)
    for line in inf:
        ll = line.strip().split("\t")
        if len(ll) <6:
            continue
        if ll[5] == '+':
            if (int(ll[1])-3) < 0 :
                print ll
                continue
            seq = genome[ll[0]][(int(ll[1])-3):(int(ll[1])+3)]
            if cutBG.has_key(seq):
                cutBG[seq]+=1
            else:
                #print seq
                pass
        elif ll[5] == '-':
            if (int(ll[2])-1-2) < 0:
                print ll
            continue
            seq = rev(genome[ll[0]][(int(ll[2])-1-2):(int(ll[2])-1+4)])
            if cutBG.has_key(seq):
                cutBG[seq]+=1
            else:
                #print seq
                pass
    inf.close()
    
    outf = open(out,'w')
    for seqtype in sorted(cutBG.keys()):
        pout = str(cutBG[seqtype]*1.0/genomeBG[seqtype])
        nout = str(cutBG[rev(seqtype)]*1.0/genomeBG[rev(seqtype)])
        #pcut = str(cut_nmer_dict_p[seqtype])
        #ncut = str(cut_nmer_dict_n[seqtype])
        #bgseq = str(allcut_nmer_dict[seqtype])
        outf.write("\t".join([seqtype,pout,nout])+"\n")
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

    optparser.add_option("-t","--tag",dest="tag",type="str",
                         help="")              
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")              
                         
    optparser.add_option("--seqbg",dest="seqbg",type="str",
                         help="6 mer sequence occurance of genome background , output of ATAC_seqbias_genome_seq.py , default = /mnt/Storage/home/huse/Project/ATAC/Data/BGmatrix/genome_seq/hg19_genome_seq_occ.txt")                                       
    optparser.add_option("--sequence",dest="sequence",type="str",default='/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("--left",dest="left",type="int",default = 3,
                         help="region length(left) for calculate seqbias , default =3  nmer = left+right")
    optparser.add_option("--right",dest="right",type="int",default = 3,
                         help="region length(right) for calculate seqbias , default =3  nmer = left+right")
        
#========minor options=============

    (options,args) = optparser.parse_args()

    tag = options.tag
    out = options.out
    seqbg = options.seqbg
    seq = options.sequence
    lflank = options.left
    rflank = options.right
    sequence = options.sequence
    if not tag or not out:
        optparser.print_help()
        sys.exit(1)
    seqbias(tag,out,sequence,lflank,rflank,seqbg)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


