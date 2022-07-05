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
    
def CGcontent(fullSEQ):
    GC=0
    for i in fullSEQ:
        if i=="G" or i=="C":
            GC +=1
#    return str(GC*1.0/(len(fullSEQ)))
    return str(GC)
        
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
        


def PEreads_feature(ll, GCrange, flank, genome2bit,biasMat):
    seqP = genome2bit[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
    seqP_shift9N = rev(genome2bit[ll[0]][(int(ll[1])-flank + 9):(int(ll[1])+flank + 9)].upper())
    seqN = rev(genome2bit[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper()) 
    seqN_shift9P = genome2bit[ll[0]][(int(ll[2])-flank - 9):(int(ll[2])+flank - 9)].upper() 
    center = (int(ll[1]) + int(ll[2]))/2
    seqGC = CGcontent(genome2bit[ll[0]][(center-GCrange/2):(center+GCrange/2)].upper())
    fraglen = int(ll[2])-int(ll[1])
    if not biasMat.has_key(seqP):#len(seqP) != flank*2 or "N" in seqP:
        biasP = "NA"
    else:
        biasP = biasMat[seqP]
    if not biasMat.has_key(seqN):#len(seqN) != flank*2 or "N" in seqN:
        biasN = "NA"
    else:
        biasN = biasMat[seqN]
    if not biasMat.has_key(seqP_shift9N):#len(seqP) != flank*2 or "N" in seqP:
        biasP_shift9N = "NA"
    else:
        biasP_shift9N = biasMat[seqP_shift9N]
    if not biasMat.has_key(seqN_shift9P):#len(seqP) != flank*2 or "N" in seqP:
        biasN_shift9P = "NA"
    else:
        biasN_shift9P = biasMat[seqN_shift9P]
    
    return ll + [biasP,biasN,seqGC,fraglen,biasP_shift9N,biasN_shift9P]

def dupFeature(inputbed, outputfilename, seq2bit, GCrange, flankLen, biasfile):
    bias = {}
    inf = open(biasfile)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        bias[ll[0]] = str(round(float(ll[2]),4))
    inf.close()
    
    genome = twobitreader.TwoBitFile(seq2bit) 
    
    inf = open(inputbed)
    outf = open(outputfilename,'w')
    
    last_reads = "NA"
    now_dup_num = 0
    read_name_idx = 0
    for line in inf:        
        if last_reads == "NA" or judge_same_read(last_reads, line):
            last_reads = line
            now_dup_num += 1
        
        else:
            # write lastreads
            read_name_idx += 1
            read_name = 'r' + str(read_name_idx)
            ll_last = last_reads.split()
            newll = PEreads_feature( ll_last[:4]+[read_name,now_dup_num], GCrange, flankLen, genome,bias)
            if not "NA" in newll:
                outf.write("\t".join(map(str,newll))+"\n")
            
            now_dup_num = 1
            last_reads = line
        
    read_name_idx += 1
    read_name = 'r' + str(read_name_idx)
    ll_last = last_reads.split()
    newll = PEreads_feature( ll_last[:4]+[read_name,now_dup_num], GCrange, flankLen, genome,bias)
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
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
    optparser.add_option("-g","--gcrange",dest="gcrange",type="int",default=100,
                         help="range for each reads to calculate GC content, default =100 (centered on reads center)")
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
    
    dupFeature(inputbed, options.outputfile, options.sequence, options.gcrange, options.flank, options.biasfile)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



