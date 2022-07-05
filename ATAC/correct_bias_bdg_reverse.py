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
import math,time
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader
import numpy
#from numpy import linalg as la

#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def rev(seq):
    revseq = ""
    for i in seq[::-1].upper():
        if i == 'A':
            r = 'T'
        elif i == 'T':
            r = 'A'
        elif i == 'C':
            r = 'G'
        elif i == 'G':
            r = 'C'
        else:
            r=i#print i
        revseq += r
    return revseq

def readBG(bgmatrix):
    BGraw = {}
    BGenc = {}
    #nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        name = ll[0]
        if ll[3] == "NA":
            BGraw[name] = "NA"
        else:
            BGraw[name] = float(ll[3])
        BGenc[name] = float(ll[4])
    inf.close()
    seqlen = len(name)
    return BGraw,BGenc,seqlen


def bias_correct_flank(inbdg,outname,biasMat,Gen,strand):

    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    
    inf = open(inbdg)
    outf_encCorrect = open(outname + ".bdg",'w')

    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        raw_sig = float(ll[3])
        if raw_sig == 0:
            outf_encCorrect.write(line)
        else:
            for pos in range(start,end):
                if strand == "+":
                    this_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
                else:
                    this_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
                if BGraw.has_key(this_seq):
                    enc_correct_sig = raw_sig / BGenc[this_seq]
                else:
                    enc_correct_sig = raw_sig
                newllenc = [chrm,pos,pos+1,enc_correct_sig]
                outf_encCorrect.write("\t".join(map(str,newllenc))+"\n")
    outf_rawCorrect.close()    
    outf_encCorrect.close()
    inf.close()

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputbdg",dest="inputbdg",type="str",
                         help="5 column summit200bp file ,sorted")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="name of output file")
    optparser.add_option("-b","--biasMat",dest="bgmatrix",type="str",default = "/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt",
                         help="bias matrix in giver N-mer, default is 8-mer estimated from Yeast naked DNA (simplex encoding)")
    optparser.add_option("-s","--strand",dest="strand",type="choice",choices=["+","-"],default = "+",
                         help="strand of bdg cut file, choose from +/-")

#========minor options=============
    optparser.add_option("--genome",dest="genome",type="str",default = "/scratch/sh8tv/Data/Genome/hg38/hg38.2bit",
                         help="genome sequence in 2bit format")


    (options,args) = optparser.parse_args()

    if not options.inputbdg:
        optparser.print_help()
        sys.exit(1)

    inputbdg = options.inputbdg
    outname = options.outname
    bgmatrix = options.bgmatrix
    gen = options.genome

    bias_correct_flank(inputbdg,outname,bgmatrix,gen,options.strand)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

