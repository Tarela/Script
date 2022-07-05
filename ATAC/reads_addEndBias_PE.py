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
    
### function

def PEreads_feature(chrm,start,end, flank, genome2bit,biasMat,TYPE):


    seqP = genome2bit[chrm][(start-flank):(start+flank)].upper()
    seqN = rev(genome2bit[chrm][(end-flank):(end+flank)].upper()) 

    if not biasMat.has_key(seqP):#len(seqP) != flank*2 or "N" in seqP:
        biasP = "NA"
    else:
        biasP = biasMat[seqP]
    if not biasMat.has_key(seqN):#len(seqN) != flank*2 or "N" in seqN:
        biasN = "NA"
    else:
        biasN = biasMat[seqN]

    if TYPE == "DNase":
        if biasP != "NA" and biasN != "NA":
            return [round(biasP,4),round(biasN,4)]
        else:
            return ["NA","NA"]
    else:
        seqP_shift9N = rev(genome2bit[chrm][(start-flank + 9):(start+flank + 9)].upper())
        seqN_shift9P = genome2bit[chrm][(end-flank - 9):(end+flank - 9)].upper() 

        if not biasMat.has_key(seqP_shift9N):#len(seqP) != flank*2 or "N" in seqP:
            biasP_shift9N = "NA"
        else:
            biasP_shift9N = biasMat[seqP_shift9N]
        if not biasMat.has_key(seqN_shift9P):#len(seqP) != flank*2 or "N" in seqP:
            biasN_shift9P = "NA"
        else:
            biasN_shift9P = biasMat[seqN_shift9P]
    
        if biasP != "NA" and biasN != "NA" and biasP_shift9N != "NA" and biasN_shift9P !="NA":
            #print biasP,biasP_shift9N,biasN,biasN_shift9P
            return [round((biasP + biasP_shift9N)/2,4),round((biasN + biasN_shift9P)/2,4)]
        else:
            return ["NA","NA"]    
#    center = (int(ll[1]) + int(ll[2]))/2
#    seqGC = CGcontent(genome2bit[ll[0]][(center-GCrange/2):(center+GCrange/2)].upper())
#    fraglen = int(ll[2])-int(ll[1])


def readsbias(inputbed, outputfilename, seq2bit, flankLen, biasfile,TYPE):
    bias = {}
    inf = open(biasfile)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        bias[ll[0]] = float(ll[2])
    inf.close()
    
    genome = twobitreader.TwoBitFile(seq2bit) 
    
    inf = open(inputbed)
    outf = open(outputfilename,'w')

    for line in inf:
#        ll = line.strip().split("\t")
        ll = line.split()
        if "." in ll[0] or "_" in ll[0]:
            continue
        if ll[0] == "chrMT":
            chrm = 'chrM'
        else:
            chrm = ll[0]

        if len(ll) < 6:
            PEtag = 1
        elif not ll[5] in ["+","-"]:
            PEtag = 1
        else:
            print "not PE bed data"
            sys.exit()
        this_bias_values = PEreads_feature(chrm,int(ll[1]),int(ll[2]),flankLen,genome,bias,TYPE)
        if "NA" in this_bias_values:
            continue
        newll = ll + this_bias_values
        outf.write("\t".join(map(str,newll))+"\n")

    outf.close()
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
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="")              
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/scratch/sh8tv/Data/Genome/hg38/hg38.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
    optparser.add_option("-b","--biasfile",dest="biasfile",type="str",default='/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt',
                         help="enc bias matrix, defualt is yeast ATAC enc 8mer bias")
    optparser.add_option("-t","--TYPE",dest="TYPE",type="str",
                         help="datatype, choose from DNase and ATAC")
 
#    optparser.add_option("-d","--dupN",dest="dupN",type="int",default=1,
#                         help="number of kept duplicate, default =1 ")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputfile = options.outputfile
    datatype  = options.TYPE
    if not inputbed or not datatype:
        optparser.print_help()
        sys.exit(1)
    if not datatype in ["ATAC", "DNase"]:
        optparser.print_help()
        sys.exit(1)

    readsbias(inputbed, options.outputfile, options.sequence, options.flank, options.biasfile, datatype)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



