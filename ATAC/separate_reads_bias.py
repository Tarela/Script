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
            r=i#print i
        revseq += r
    return revseq   


def PEreads_feature(ll, flank, genome2bit, biasMat, datatype):
    seqP = genome2bit[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
    seqN = rev(genome2bit[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper()) 
    if not biasMat.has_key(seqP):#len(seqP) != flank*2 or "N" in seqP:
        biasP = "NA"
    else:
        biasP = biasMat[seqP]
    if not biasMat.has_key(seqN):#len(seqN) != flank*2 or "N" in seqN:
        biasN = "NA"
    else:
        biasN = biasMat[seqN]

    if datatype=="ATAC":
        seqP_shift9N = rev(genome2bit[ll[0]][(int(ll[1])-flank + 9):(int(ll[1])+flank + 9)].upper())
        seqN_shift9P = genome2bit[ll[0]][(int(ll[2])-flank - 9):(int(ll[2])+flank - 9)].upper() 

        if not biasMat.has_key(seqP_shift9N):#len(seqP) != flank*2 or "N" in seqP:
            biasP_shift9N = "NA"
        else:
            biasP_shift9N = biasMat[seqP_shift9N]
        if not biasMat.has_key(seqN_shift9P):#len(seqP) != flank*2 or "N" in seqP:
            biasN_shift9P = "NA"
        else:
            biasN_shift9P = biasMat[seqN_shift9P]

        if biasP == "NA" or biasP_shift9N == "NA":
            biasLeft = "NA"
        else:
            biasLeft = ( biasP + biasP_shift9N ) /2
        if biasN == "NA" or biasN_shift9P == "NA":
            biasRight = "NA"
        else:
            biasRight = ( biasN + biasN_shift9P ) /2

        return biasLeft,biasRight

    else:
        return biasP,biasN

def dupFeature(usename, seq2bit, datatype):

    if datatype == "ATAC":
        biasfile = "/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt"
        cutoff = 3.5
    elif datatype == "DNase":
        biasfile = "/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedIMR90_DNase_Enc8mer.txt"
        cutoff = -2.2
    else:
        print datatype
        sys.exit()

    bias = {}
    inf = open(biasfile)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        bias[ll[0]] = round(float(ll[2]),4)
    inf.close()
    flankLen = len(ll[0])/2
    genome = twobitreader.TwoBitFile(seq2bit) 
    
    inf = open(usename + "_uniq.bed")
    outf1 = open(usename + "_g1.bed",'w')
    outf2 = open(usename + "_g2.bed",'w')
    outf3 = open(usename + "_g3.bed",'w')
    outf4 = open(usename + "_g4.bed",'w')
    
    for line in inf:        
        ll = line.split()
        if ll[0] == "chrM":
            continue

        left,right = PEreads_feature( ll, flankLen, genome, bias, datatype)

        if left == "NA" or right == "NA" :   
            continue

        newll = [ll[0],( int(ll[1]) + int(ll[2]) )/2, ( int(ll[1]) + int(ll[2]) )/2 + 1, ".",".","+"]
        newline =  "\t".join( map(str,newll) ) + "\n" 

        if left > cutoff:
            if right > cutoff:
                outf1.write(newline)
            else:
                outf4.write(newline)
        else:
            if right > cutoff:
                outf2.write(newline)
            else:
                outf3.write(newline)

    outf1.close()
    outf2.close()
    outf3.close()
    outf4.close()

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
          
    optparser.add_option("-n","--name",dest="usename",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/scratch/sh8tv/Data/Genome/hg38/hg38.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-t","--datatype",dest="datatype",type="str",
                         help="choose from DNase and ATAC")              
 
#    optparser.add_option("-d","--dupN",dest="dupN",type="int",default=1,
#                         help="number of kept duplicate, default =1 ")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    usename = options.usename
    if not usename:
        optparser.print_help()
        sys.exit(1)
    dupFeature(usename, options.sequence, options.datatype)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



