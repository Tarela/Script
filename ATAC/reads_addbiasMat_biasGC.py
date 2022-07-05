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
    

#    return str(GC)
def seq2GC(inputseq):
    GC = 0
    for i in inputseq:
        if i.upper()=="G" or i.upper()=="C":
            GC +=1
    return GC*1.0/len(inputseq)
       
def judge_same_read(readA,readB):
    llA = readA.split()
    llB = readB.split()
    if llA[0] == llB[0] and llA[1] == llB[1] and llA[2] == llB[2] :
        return 1
    else:
        return 0


### function
def collectBias(infile,outfile,biasMat):
    inf = open(infile)
    outf = open(outfile,'w')
    for line in inf:
        ll = line.strip().split("\t")
        if len(ll) == 8:
            if not ll[5] in ["+","-"]:
                print "SE strand bug",outfile,len(ll),ll
                continue
            if len(ll[6].strip()) != flank*2 or "N" in ll[6].strip():
                print "seqlength bug",outfile,len(ll[6].strip()),ll
                continue
            newll = ll[:6] +  [biasMat[ll[6]],ll[7]]
            outf.write("\t".join(newll)+"\n")
        elif len(ll) == 9:
            if len(ll[5].strip()) != flank*2 or len(ll[6].strip()) != flank*2 or "N" in ll[5].strip() or "N" in ll[6].strip():
                print "seqlength bug",outfile,len(ll[5].strip()),len(ll[6].strip()),ll
                continue
            newll = ll[:5] +  [biasMat[ll[5]], biasMat[ll[6]] ,ll[7],ll[8]]
            outf.write("\t".join(newll)+"\n")
        else:
            print "file colNum bug",outfile,len(ll),ll
        
def CBvalue(V1,V2):
    if V1 =="NA" or V2=="NA":
        return "NA"
    else:
        #print V1,V2
        return (V1+V2)/2

def seq2GC(inputseq):
    GC = 0
    for i in inputseq:
        if i.upper()=="G" or i.upper()=="C":
            GC +=1
    return GC*1.0/len(inputseq)


def SEreads_feature(inline, GCrange, flank, genome2bit,biasMat):
    ll = inline.split()
    if ll[5] == "+":
        seq = genome2bit[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
    #seqP_shift9N = rev(genome2bit[ll[0]][(int(ll[1])-flank + 9):(int(ll[1])+flank + 9)].upper())
        readseq = genome2bit[ll[0]][(int(ll[1])):(int(ll[1])+GCrange)].upper()
    else:   
        seq = rev(genome2bit[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper()) 
        readseq = rev(genome2bit[ll[0]][(int(ll[2])-GCrange):(int(ll[2]))].upper()) 
    #seqN_shift9P = genome2bit[ll[0]][(int(ll[2])-flank - 9):(int(ll[2])+flank - 9)].upper() 
    #center = (int(ll[1]) + int(ll[2]))/2
    #seqGC = seq2GC(genome2bit[ll[0]][(center-GCrange/2):(center+GCrange/2)].upper())
    #fraglen = int(ll[2])-int(ll[1])
    if not biasMat.has_key(seq):#len(seqP) != flank*2 or "N" in seqP:
        bias = "NA"
        biasGC = "NA"
        readsGC = "NA"
    else:
        bias = biasMat[seq]
        biasGC = seq2GC(seq)
        readsGC = seq2GC(readseq)
#    if not biasMat.has_key(seqN):#len(seqN) != flank*2 or "N" in seqN:
#        biasN = "NA"
#    else:
#        biasN = biasMat[seqN]
#    if not biasMat.has_key(seqP_shift9N):#len(seqP) != flank*2 or "N" in seqP:
#        biasP_shift9N = "NA"
#    else:
#        biasP_shift9N = biasMat[seqP_shift9N]
#    if not biasMat.has_key(seqN_shift9P):#len(seqP) != flank*2 or "N" in seqP:
#        biasN_shift9P = "NA"
#    else:
#        biasN_shift9P = biasMat[seqN_shift9P]
    return ll + [bias,biasGC,readsGC]#biasN,CBvalue(biasP,biasP_shift9N),CBvalue(biasN,biasN_shift9P),seqGC]

def dupFeature(inputbed, outputfilename, seq2bit, GCrange, biasfile):
    bias = {}
    inf = open(biasfile)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        bias[ll[0]] = (round(float(ll[2]),4))
        flankLen= len(ll[0])
    inf.close()
    
    genome = twobitreader.TwoBitFile(seq2bit) 
    
    inf = open(inputbed)
    outf = open(outputfilename,'w')
    
    for line in inf:        
        newll = SEreads_feature( line, GCrange, flankLen, genome,bias)
        if not "NA" in newll:
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
    optparser.add_option("-g","--gcrange",dest="gcrange",type="int",default=50,
                         help="range for each reads to calculate GC content, default =50 (centered on reads center)")
    optparser.add_option("-b","--biasfile",dest="biasfile",type="str",default='/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt',
                         help="enc bias matrix, defualt is yeast ATAC enc 8mer bias")
 
#    optparser.add_option("-d","--dupN",dest="dupN",type="int",default=1,
#                         help="number of kept duplicate, default =1 ")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputfile = options.outputfile
    if not inputbed:
        optparser.print_help()
        sys.exit(1)
    
    dupFeature(inputbed, options.outputfile, options.sequence, options.gcrange, options.biasfile)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



