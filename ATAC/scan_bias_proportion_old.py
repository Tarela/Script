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
            r=i#print i
        revseq += r
    return revseq

def readBG(bgmatrix):
    BG = {}
    #nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        ll = line.split()
        name = ll[0]
        seqlen = len(name)
        BG[name] = float(ll[2])
        #nBG[name] = float(ll[2])
    inf.close()
    return BG,seqlen


def sitepro_scan(peak,outname,biasMat,Cspan,Gen,offset):

#    t = time.time()
    genome = twobitreader.TwoBitFile(Gen)
    
    BG,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    
    inf = open(peak)
#    w_plus_H=BigWigFile(open(w_plus, 'rb'))
#    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    outf_rawPlus = open(outname + "_rawPlus.bdg",'w')
    outf_rawMinus = open(outname + "_rawMinus.bdg",'w')
    outf_cbPlus = open(outname + "_cbPlus.bdg",'w')
    outf_cbMinus = open(outname + "_cbMinus.bdg",'w')
    
    outf_rawPlus_prop = open(outname + "_rawPlusProp.bdg",'w')
    outf_rawMinus_prop = open(outname + "_rawMinusProp.bdg",'w')
    outf_cbPlus_prop = open(outname + "_cbPlusProp.bdg",'w')
    outf_cbMinus_prop = open(outname + "_cbMinusProp.bdg",'w')
    
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
    
        plus_single_bias_vector = []
        plus_cb_bias_vector = []
        minus_single_bias_vector = []
        minus_cb_bias_vector = []
    
        for pos in range(start-Cspan,end+Cspan):
            plus_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
            plus_reverse_seq = rev(genome[chrm][(pos+offset-flank):(pos+offset+flank)].upper())
            minus_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
            minus_reverse_seq = genome[chrm][(pos-offset+1-flank):(pos-offset+1+flank)].upper()

            ### calculate bias value for each bases
            
            if len(plus_seq) == Nmer and not "N" in plus_seq:
                plus_bias = BG[plus_seq]
            else:
                plus_bias = -10
                
            if len(plus_reverse_seq) == Nmer and not "N" in plus_reverse_seq:
                plus_reverse_bias = BG[plus_reverse_seq]
            else:
                plus_reverse_bias = -10
            
            if len(minus_seq) == Nmer and not "N" in minus_seq:
                minus_bias = BG[minus_seq]
            else:
                minus_bias = -10

            if len(minus_reverse_seq) == Nmer and not "N" in minus_reverse_seq:
                minus_reverse_bias = BG[minus_reverse_seq]
            else:
                minus_reverse_bias = -10

            if plus_bias != "NA" and plus_reverse_bias != "NA":
                plus_cb_bias = (plus_bias + plus_reverse_bias ) / 2
            else:
                plus_cb_bias = -10
            
            if minus_bias != "NA" and minus_reverse_bias != "NA":
                minus_cb_bias = (minus_bias + minus_reverse_bias ) / 2
            else:
                minus_cb_bias = -10

            plus_single_bias_vector.append(plus_bias)
            plus_cb_bias_vector.append(plus_cb_bias)
            minus_single_bias_vector.append(minus_bias)
            minus_cb_bias_vector.append(minus_cb_bias)

        #if "NA" in plus_single_bias_vector or "NA" in plus_cb_bias_vector or "NA" in minus_single_bias_vector or "NA" in minus_cb_bias_vector:
        #    print ll
        #    continue        
        #### construct bias vector and transform to linear
        Plus_Single_Bias_log = numpy.array(plus_single_bias_vector)
        Plus_Combine_Bias_log = numpy.array(plus_cb_bias_vector)
        Minus_Single_Bias_log = numpy.array(minus_single_bias_vector)
        Minus_Combine_Bias_log = numpy.array(minus_cb_bias_vector)

        Plus_Single_Bias_linear = pow(numpy.e,Plus_Single_Bias_log)
        Plus_Combine_Bias_linear = pow(numpy.e,Plus_Combine_Bias_log)
        Minus_Single_Bias_linear = pow(numpy.e,Minus_Single_Bias_log)
        Minus_Combine_Bias_linear = pow(numpy.e,Minus_Combine_Bias_log)

        roundN = 4

        #### assign bias to bp and proportion
        for outpos in range(Cspan,(end-start+Cspan)):
            this_plus_single = round(Plus_Single_Bias_log[outpos],roundN)
            this_plus_cb = round(Plus_Combine_Bias_log[outpos],roundN)
            this_minus_single = round(Minus_Single_Bias_log[outpos],roundN)
            this_minus_cb = round(Minus_Combine_Bias_log[outpos],roundN)

            this_plus_single_prop = round(Plus_Single_Bias_linear[outpos]/sum(Plus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            this_plus_cb_prop = round(Plus_Combine_Bias_linear[outpos]/sum(Plus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            this_minus_single_prop = round(Minus_Single_Bias_linear[outpos]/sum(Minus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            this_minus_cb_prop = round(Minus_Combine_Bias_linear[outpos]/sum(Minus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)

            out_chrm = chrm
            out_start = start + outpos - Cspan
            out_end = out_start+1
            
            outf_rawPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single]))+"\n")
            outf_rawMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single]))+"\n")
            outf_cbPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb]))+"\n")
            outf_cbMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb]))+"\n")
            
            outf_rawPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single_prop]))+"\n")
            outf_rawMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single_prop]))+"\n")
            outf_cbPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb_prop]))+"\n")
            outf_cbMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb_prop]))+"\n")

    outf_rawPlus.close()
    outf_rawMinus.close()
    outf_cbPlus.close()
    outf_cbMinus.close()
    outf_rawPlus_prop.close()
    outf_rawMinus_prop.close()
    outf_cbPlus_prop.close()
    outf_cbMinus_prop.close()

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
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="5 column summit200bp file ,sorted")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="name of output file")
    optparser.add_option("-b","--biasMat",dest="bgmatrix",type="str",default = "/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt",
                         help="bias matrix in giver N-mer, default is 8-mer estimated from Yeast naked DNA (simplex encoding)")

#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 25,
                         help="region for get total signal in single bp, default = 25 means +-25bp(total 50bp) signal as total for each bp")
    optparser.add_option("--genome",dest="genome",type="str",default = "/scratch/sh8tv/Data/Genome/hg38/hg38.2bit",
                         help="genome sequence in 2bit format")
    optparser.add_option("--offset",dest="offset",type="int",default = 9,
                         help="offset related, distance of pair of +/- related cut,default = 9")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    outname = options.outname
    bgmatrix = options.bgmatrix
    gen = options.genome
    Cspan = options.Cspan
    offset = options.offset
    sitepro_scan(interval,outname,bgmatrix,Cspan,gen,offset)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

