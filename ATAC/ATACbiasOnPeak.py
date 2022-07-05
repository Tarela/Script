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
            r=i#print i
        revseq += r
    return revseq
    
def CGcontent(fullSEQ):
    GC=0
    for i in fullSEQ:
        if i=="G" or i=="C":
            GC +=1
    return str(GC*1.0/(len(fullSEQ)))
#    return str(GC)
        
### function        

def addavebias(inputbed, outputfilename, seq2bit, biasfile, plusCutbw, minusCutbw, biasmode):
    bias = {}
    inf = open(biasfile)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        bias[ll[0]] = (round(float(ll[2]),4))
    inf.close()
    flank = len(ll[0])/2
    
    genome2bit = twobitreader.TwoBitFile(seq2bit) 
    w_plus_H=BigWigFile(open(plusCutbw, 'rb'))
    w_minus_H=BigWigFile(open(minusCutbw, 'rb'))

    inf = open(inputbed)
    outf = open(outputfilename,'w')
    
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        #fullseq = genome2bit[chrm][start:end].upper()
        #peakGC = CGcontent(fullseq)
        Pcuts = list(w_plus_H.summarize(chrm,start,end,end-start).sum_data)
        Ncuts = list(w_minus_H.summarize(chrm,start,end,end-start).sum_data)
        total_cuts = 0
        total_bias = 0
        for position in range(start,end):
            cut_p = Pcuts[position-start]
            cut_n = Ncuts[position-start]
            seqP = genome2bit[chrm][(position-flank):(position+flank)].upper()
            seqP_shift9N = rev(genome2bit[chrm][(position-flank + 9):(position+flank + 9)].upper())
            seqN = rev(genome2bit[chrm][(position+1-flank):(position+1+flank)].upper()) 
            seqN_shift9P = genome2bit[chrm][(position+1-flank - 9):(position+1+flank - 9)].upper()
            if not bias.has_key(seqP) or not bias.has_key(seqP_shift9N): 
                print ll,position,"+",seqP,seqP_shift9N
                bias_p = 0
            else:
                if biasmode == "FR":
                    bias_p = int(cut_p)*((bias[seqP]+bias[seqP_shift9N])/2)
                else:
                    bias_p = int(cut_p)*(bias[seqP])
            if not bias.has_key(seqN) or not bias.has_key(seqN_shift9P): 
                print ll,position,"-",seqN,seqN_shift9P
                bias_n = 0
            else:
                if biasmode == "FR":
                    bias_n = int(cut_n)*((bias[seqN]+bias[seqN_shift9P])/2)
                else:
                    bias_n = int(cut_n)*(bias[seqN])
            total_cuts += cut_p + cut_n
            total_bias += bias_p + bias_n
        if total_cuts == 0:
            print ll,"total_cuts = 0"
            ave_bias = 0
        else:
            ave_bias = total_bias/total_cuts
        newll = ll + [ave_bias, total_cuts]

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
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")              
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/scratch/sh8tv/Data/Genome/hg38/hg38.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-b","--biasfile",dest="biasfile",type="str",default='/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt',
                         help="enc bias matrix, defualt is yeast ATAC enc 8mer bias")
    optparser.add_option("-p","--pluscutbw",dest="pluscutbw",type="str",
                         help="bigwig file of plus strand cuts")              
    optparser.add_option("-n","--minuscutbw",dest="minuscutbw",type="str",
                         help="bigwig file of minus strand cuts")              
    optparser.add_option("-m","--biasmode",dest="biasmode",type="str",default="F",
                         help="bias mode, choose from F or FR, represent forward only or forwardXreverse")
#    optparser.add_option("-d","--dupN",dest="dupN",type="int",default=1,
#                         help="number of kept duplicate, default =1 ")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputfile = options.outputfile
    pluscut = options.pluscutbw
    minuscut = options.minuscutbw
    biasmode = options.biasmode
    if not inputbed or not outputfile or not pluscut or not minuscut:
        optparser.print_help()
        sys.exit(1)
    if not biasmode in ["F","FR"]:
        optparser.print_help()
        sys.exit(1)

    addavebias(inputbed, outputfile, options.sequence, options.biasfile, pluscut, minuscut, biasmode)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



