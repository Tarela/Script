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
            BGraw[name] = 0
        else:
            BGraw[name] = pow(numpy.e,float(ll[1]))
        BGenc[name] = pow(numpy.e,float(ll[2]))
    inf.close()
    seqlen = len(name)
    return BGraw,BGenc,seqlen


def bias_scan_fxr(peak,outname,biasMat,Gen):

#    t = time.time()
    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    
    inf = open(peak)
#    w_plus_H=BigWigFile(open(w_plus, 'rb'))
#    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf1 = open(outname+"_rawbias.bed",'w')
    outf2 = open(outname+"_encbias.bed",'w')
#    outf3 = open(outname+"_fxrbias.bed",'w')
    
    for line in inf:
        #ll = line.split()
        #chrm = ll[0]
        #start = int(ll[1])
        #end = int(ll[2])
        ll = line.split()
        if ll[0] == "chrM" or "_" in ll[0]:
            continue
        motiflen = int(ll[2]) - int(ll[1])
        #upstream_ext = fulllen/2 - motiflen/2
        if ll[5] == "+":
            center = int(ll[1]) + motiflen/2
        else:
            center = int(ll[2]) - motiflen/2 
        #center = int((int(ll[2]) + int(ll[1]))/2)
        chrm = ll[0]
        start = center-100#int(ll[1])-100#upstream_ext
        end = center+100#start + 200#fulllen

        plus_single_bias_raw_vector = []
        plus_single_bias_enc_vector = []
#        plus_single_bias_fxr_vector = []
        minus_single_bias_raw_vector = []
        minus_single_bias_enc_vector = []
#        minus_single_bias_fxr_vector = []
    
        #for pos in range(start-Cspan,end+Cspan):
        for pos in range(start,end):
            plus_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
#            plus_reverse_seq = rev(genome[chrm][(pos+offset-flank):(pos+offset+flank)].upper())
            minus_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
#            minus_reverse_seq = genome[chrm][(pos-offset+1-flank):(pos-offset+1+flank)].upper()

            ### calculate bias value for each bases
            
            if len(plus_seq) == Nmer and not "N" in plus_seq:
                plus_bias_raw = BGraw[plus_seq]
                plus_bias_enc = BGenc[plus_seq]
            else:
                plus_bias_raw = 0
                plus_bias_enc = 0

#            if len(plus_reverse_seq) == Nmer and not "N" in plus_reverse_seq:
#                plus_reverse_bias = BGenc[plus_reverse_seq]
#            else:
#                plus_reverse_bias = 0
            
            if len(minus_seq) == Nmer and not "N" in minus_seq:
                minus_bias_raw = BGraw[minus_seq]
                minus_bias_enc = BGenc[minus_seq]
            else:
                minus_bias_raw = 0
                minus_bias_enc = 0

 #           if len(minus_reverse_seq) == Nmer and not "N" in minus_reverse_seq:
 #               minus_reverse_bias = BGenc[minus_reverse_seq]
 #           else:
 #               minus_reverse_bias = 0

 #           plus_cb_bias = numpy.sqrt(plus_bias_enc * plus_reverse_bias ) 
 #           minus_cb_bias = numpy.sqrt(minus_bias_enc * minus_reverse_bias)

            plus_single_bias_raw_vector.append(round(plus_bias_raw,4))
            plus_single_bias_enc_vector.append(round(plus_bias_enc,4))
#            plus_single_bias_fxr_vector.append(round(plus_cb_bias,4))
            #plus_cb_bias_vector.append(plus_cb_bias)
            minus_single_bias_raw_vector.append(round(minus_bias_raw,4))
            minus_single_bias_enc_vector.append(round(minus_bias_enc,4))
#            minus_single_bias_fxr_vector.append(round(minus_cb_bias,4))
            #minus_cb_bias_vector.append(minus_cb_bias)

        #### construct bias vector and transform to linear

        if ll[5] == "+":
            Plus_Single_rawBias = list(plus_single_bias_raw_vector)
            Plus_Single_encBias = list(plus_single_bias_enc_vector)
#            Plus_Single_fxrBias = list(plus_single_bias_fxr_vector)
            Minus_Single_rawBias = list(minus_single_bias_raw_vector)
            Minus_Single_encBias = list(minus_single_bias_enc_vector)
#            Minus_Single_fxrBias = list(minus_single_bias_fxr_vector)
        else:
            Plus_Single_rawBias = list(minus_single_bias_raw_vector)[::-1]
            Plus_Single_encBias = list(minus_single_bias_enc_vector)[::-1]
#            Plus_Single_fxrBias = list(minus_single_bias_fxr_vector)[::-1]
            Minus_Single_rawBias = list(plus_single_bias_raw_vector)[::-1]
            Minus_Single_encBias = list(plus_single_bias_enc_vector)[::-1]
#            Minus_Single_fxrBias = list(plus_single_bias_fxr_vector)[::-1]

        newll1 = ll + Plus_Single_rawBias + Minus_Single_rawBias
        outf1.write("\t".join(map(str,newll1))+"\n")
        newll2 = ll + Plus_Single_encBias + Minus_Single_encBias
        outf2.write("\t".join(map(str,newll2))+"\n")
#        newll3 = ll + Plus_Single_fxrBias + Minus_Single_fxrBias
#        outf3.write("\t".join(map(str,newll3))+"\n")

    inf.close()
    outf1.close()
    outf2.close()
    #outf3.close()




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
                         help="motif 6 column ,sorted")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="name of output file")
    optparser.add_option("-b","--biasMat",dest="bgmatrix",type="str",default = "/nv/vol190/zanglab/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedIMR90_DNase_Enc8mer.txt",
                         help="bias matrix in giver N-mer, default is 8-mer estimated from Yeast naked DNA (simplex encoding)")
    #optparser.add_option("-t","--biastype",dest="biastype",type="choice",choices=["sob","flank","fxr"],default = "flank",
    #                     help="bias type, choices from flank, sob and fxr")

#========minor options=============
    #optparser.add_option("--Cspan",dest="Cspan",type="int",default = 100,
    #                     help="region for get total signal in single bp, default = 100 means +-100bp(total 200bp) signal as total for each bp")
    optparser.add_option("--genome",dest="genome",type="str",default = "/nv/vol190/zanglab/sh8tv/Data/Genome/hg38/hg38.2bit",
                         help="genome sequence in 2bit format")
    #optparser.add_option("--offset",dest="offset",type="int",default = 9,
    #                     help="offset related, distance of pair of +/- related cut,default = 9")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    outname = options.outname
    bgmatrix = options.bgmatrix
    gen = options.genome
    #Cspan = options.Cspan
    #offset = options.offset
    #biastype = options.biastype
    #if biastype == "flank":
      #  bias_scan_flank(interval,outname,bgmatrix,Cspan,gen)
    #if biastype == "sob":
     #   bias_scan_sob(interval,outname,bgmatrix,Cspan,gen)
    #if biastype == "fxr":
    bias_scan_fxr(interval,outname,bgmatrix,gen)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

