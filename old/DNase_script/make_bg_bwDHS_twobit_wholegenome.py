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
#from cistrome import regions as c
import twobitreader
from bx.bbi.bigwig_file import BigWigFile
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
def count_cut_nmers( cut_plus,cut_minus,lflank, rflank ,single_nmer_cutoff,sequence):
    """
    count the number of cuts associated with each nmer in sequence covered by X.
    offset is the position of the cut to be associated with each nmer.
    if offset = 0 the first base of the tag is lined up with the nmer start
    """
    #w_plus_H=BigWigFile(open(w_plus, 'rb'))
    #w_minus_H=BigWigFile(open(w_minus, 'rb'))

    genome = twobitreader.TwoBitFile(sequence)
    # keep count of the number of occurrences of each n-mer
    
    seq_nmer_dict   = make_nmer_dict(rflank+lflank)
    for chrom in genome.keys():
        wholeSeq = genome[chrom][:]
        RVwholeSeq = rev(wholeSeq)
        for i in range(len(wholeSeq)-lflank-rflank +1 ):
            seq_nmer_dict[wholeSeq[i:(i+lflank+rflank)]] += 1
            seq_nmer_dict[RVwholeSeq[i:(i+lflank+rflank)]] += 1
    #cut_nmer_dict    = {}

    cut_nmer_dict = make_nmer_dict(rflank+lflank)
    inf = open(cut_plus)
    for line in inf:
        ll = line.strip().split("\t")
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        value = float(ll[3])
        if value == 0:
            continue

        for i in range(start,end):   
            try:     
                seq = genome[chrm][(i-lflank):(i+rflank)]
            except:
                print 'plus',chrm,start,end,value
                continue
            if not cut_nmer_dict.has_key(seq):
                continue	
            if value <= single_nmer_cutoff:
                cut_nmer_dict[seq] += value
            else:
                seq_nmer_dict[seq] -= 1
    inf.close()
    inf = open(cut_minus)
    for line in inf:
        ll = line.strip().split("\t")
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        value = float(ll[3])
        if value == 0:
            continue
        for i in range(start,end):  
            try :      
                seq = rev(genome[chrm][(i-rflank+1):(i+lflank+1)])
            except : 
                print 'minus',chrm,start,end,value
                continue
            if not cut_nmer_dict.has_key(seq):
                continue
            if value <= single_nmer_cutoff:
                cut_nmer_dict[seq] += value
            else:
                seq_nmer_dict[seq] -= 1
            
    inf.close()

    return seq_nmer_dict, cut_nmer_dict



def build_nmer_model(cut_plus,cut_minus,out,single_nmer_cut,sequence,lflank,rflank):
    """
    count cuts within each nmer sequence 
    """
    # |0.1.2.3.4 lflank=0 rflank=5
    #  0|1.2.3.4 lflank=1 rflank=4
    nmer = lflank + rflank

    # count  
    seq_nmer_dict, cut_nmer_dict = count_cut_nmers( cut_plus,cut_minus, lflank,rflank ,single_nmer_cut,sequence)
    # -1 : no seq , 0 : seq no cut
    dall = make_nmer_dict(nmer)
    #print len(allcut_nmer_dict),len(cut_nmer_dict_p),len(cut_nmer_dict_n)
    # normalize
    for elem in seq_nmer_dict:

        dall[elem]  =  (1.0*cut_nmer_dict[elem]) / (1.0*seq_nmer_dict[elem]) 
    outf = open(out,'w')
    for seqtype in sorted(dall.keys()):
        pout = str(dall[seqtype])
        nout = str(dall[rev(seqtype)])
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
    optparser.add_option("-o","--output",dest="out",type="str",
                         help="")
    optparser.add_option("-p","--cutplus",dest="cutp",type="str",
                         help="cleavage bedGraph file for DNase positive cut")
    optparser.add_option("-n","--cutminus",dest="cutm",type="str",
                         help="cleavage bedGraph file for DNase negative cut")
    
    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/Scerevisiae.2bit',
                         help="twobit sequence file , default : /home/sh430/Data/Genome/Scerevisiae.2bit")
    optparser.add_option("--sc",dest="single_nmer_cut",type="int",default = 20,
                         help="cutoff for cleavage count at single loci , if there's more than N cut at a single loci , we think its caused by amplification ,and skip it")
    optparser.add_option("--left",dest="left",type="int",default = 3,
                         help="region length(left) for calculate seqbias , default =3  nmer = left+right")
    optparser.add_option("--right",dest="right",type="int",default = 3,
                         help="region length(right) for calculate seqbias , default =3  nmer = left+right")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    out = options.out
    cut_plus = options.cutp
    cut_minus = options.cutm
    genome = options.genome
    single_nmer_cut = options.single_nmer_cut

    if not cut_plus or not cut_minus:
        optparser.print_help()
        sys.exit(1)
    
    build_nmer_model(cut_plus,cut_minus,out,single_nmer_cut,genome,options.left,options.right)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



