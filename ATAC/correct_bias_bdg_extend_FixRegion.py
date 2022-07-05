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
        if ll[1] == "NA":
            BGraw[name] = "NA"
        else:
            BGraw[name] = float(ll[3])
        BGenc[name] = float(ll[4])
    inf.close()
    seqlen = len(name)
    return BGraw,BGenc,seqlen


def bias_correct_flank(biasMat,Gen):

    ext = 50
    outname = "GM12878_ATAC"
    inbdg_plus = "/scratch/sh8tv/Project/scATAC/Result/diffmer_var_comparison/obs_pred_reads_comparison/cleavage/GM12878_SEATACr1_PlusCuts.bdg"
    inbdg_minus = "/scratch/sh8tv/Project/scATAC/Result/diffmer_var_comparison/obs_pred_reads_comparison/cleavage/GM12878_SEATACr1_MinusCuts.bdg"
    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    
    outf_raw = open(outname + "_rawsig_ext%s.bdg"%ext,'w')
    outf_correct = open(outname + "_encFlank_ext%s.bdg"%ext,'w')

    infp = open(inbdg_plus)
    infn = open(inbdg_minus)

    orgstart = 141305021
    orgend = 141309021
    orgchrm = "chr4"
    region_len = (orgend - orgstart)
    rawlist = numpy.array([0.0]*region_len)
    correctlist = numpy.array([0.0]*region_len)
    for line in infp:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        
        if chrm == orgchrm and start >= orgstart and start <= orgend:
            pass
        else:
            continue
        raw_sig = float(ll[3])
        for pos in range(start,end):
            this_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
            if BGenc.has_key(this_seq):
                enc_correct_sig = raw_sig * BGenc[this_seq]
            else:
                enc_correct_sig = raw_sig
            idx_start = pos - orgstart
            idx_end = min(idx_start+ext,region_len)
            rawlist[idx_start:idx_end] += raw_sig
            correctlist[idx_start:idx_end] += enc_correct_sig
    infp.seek(0)

    for line in infn:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        
        if chrm == orgchrm and start >= orgstart and start <= orgend:
            pass
        else:
            continue
        raw_sig = float(ll[3])
        for pos in range(start,end):
            this_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
            if BGenc.has_key(this_seq):
                enc_correct_sig = raw_sig * BGenc[this_seq]
            else:
                enc_correct_sig = raw_sig
            idx_start = max(0,pos - ext - orgstart)
            idx_end = pos - orgstart
            rawlist[idx_start:idx_end] += raw_sig
            correctlist[idx_start:idx_end] += enc_correct_sig
    infn.seek(0)

    for i in range(len(rawlist)):
        newll = [orgchrm, orgstart + i,orgstart+i+1,rawlist[i]]
        outf_raw.write("\t".join(map(str,newll))+"\n")

    for i in range(len(correctlist)):
        newll = [orgchrm, orgstart + i,orgstart+i+1,correctlist[i]]
        outf_correct.write("\t".join(map(str,newll))+"\n")
   


    orgstart = 117480877
    orgend = 117484877
    orgchrm = "chr6"
    region_len = (orgend - orgstart)
    rawlist = numpy.array([0.0]*region_len)
    correctlist = numpy.array([0.0]*region_len)
    for line in infp:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        
        if chrm == orgchrm and start >= orgstart and start <= orgend:
            pass
        else:
            continue
        raw_sig = float(ll[3])
        for pos in range(start,end):
            this_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
            if BGenc.has_key(this_seq):
                enc_correct_sig = raw_sig * BGenc[this_seq]
            else:
                enc_correct_sig = raw_sig
            idx_start = pos - orgstart
            idx_end = min(idx_start+ext,region_len)
            rawlist[idx_start:idx_end] += raw_sig
            correctlist[idx_start:idx_end] += enc_correct_sig
    infp.close()

    for line in infn:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        
        if chrm == orgchrm and start >= orgstart and start <= orgend:
            pass
        else:
            continue
        raw_sig = float(ll[3])
        for pos in range(start,end):
            this_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
            if BGenc.has_key(this_seq):
                enc_correct_sig = raw_sig * BGenc[this_seq]
            else:
                enc_correct_sig = raw_sig
            idx_start = max(0,pos - ext - orgstart)
            idx_end = pos - orgstart
            rawlist[idx_start:idx_end] += raw_sig
            correctlist[idx_start:idx_end] += enc_correct_sig
    infn.close()

    for i in range(len(rawlist)):
        newll = [orgchrm, orgstart + i,orgstart+i+1,rawlist[i]]
        outf_raw.write("\t".join(map(str,newll))+"\n")

    for i in range(len(correctlist)):
        newll = [orgchrm, orgstart + i,orgstart+i+1,correctlist[i]]
        outf_correct.write("\t".join(map(str,newll))+"\n")

    outf_raw.close()
    outf_correct.close()    



# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-b","--biasMat",dest="bgmatrix",type="str",default = "/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt",
                         help="bias matrix in giver N-mer, default is 8-mer estimated from Yeast naked DNA (simplex encoding)")

#========minor options=============
    optparser.add_option("--genome",dest="genome",type="str",default = "/scratch/sh8tv/Data/Genome/hg38/hg38.2bit",
                         help="genome sequence in 2bit format")


    (options,args) = optparser.parse_args()

    if not options.bgmatrix:
        optparser.print_help()
        sys.exit(1)

    bgmatrix = options.bgmatrix
    gen = options.genome

    bias_correct_flank(bgmatrix,gen)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

