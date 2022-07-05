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
    if llA[0] == llB[0] and llA[1] == llB[1] and llA[2] == llB[2] and llA[5] == llB[5]:
        return 1
    else:
        return 0

def SEreads_feature(ll, GCrange, flank, genome2bit):
    Strand = ll[5]
    if Strand == "+":
        seq = genome2bit[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
    else:
        seq = rev(genome2bit[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper()) 
    center = (int(ll[1]) + int(ll[2]))/2
    seqGC = CGcontent(genome2bit[ll[0]][(center-GCrange/2):(center+GCrange/2)].upper())
    #fraglen = int(ll[2])-int(ll[1])
    return ll + [seq,seqGC]

def dupFeature(inputbed, outputfilename, seq2bit, GCrange, flankLen):
    
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
            newll = SEreads_feature( ll_last[:3]+[read_name,now_dup_num,ll_last[5]], GCrange, flankLen, genome)
            outf.write("\t".join(map(str,newll))+"\n")
            
            now_dup_num = 1
            last_reads = line
        
    read_name_idx += 1
    read_name = 'r' + str(read_name_idx)
    ll_last = last_reads.split()
    newll = SEreads_feature( ll_last[:3]+[read_name,now_dup_num,ll_last[5]], GCrange, flankLen, genome)
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
    optparser.add_option("-o","--out",dest="outputfile",type="str",
                         help="")
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/nv/vol190/zanglab/sh8tv/Data/Genome/hg38/hg38.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
    optparser.add_option("-g","--gcrange",dest="gcrange",type="int",default=100,
                         help="range for each reads to calculate GC content, default =100 (centered on reads center)")
 
#    optparser.add_option("-d","--dupN",dest="dupN",type="int",default=1,
#                         help="number of kept duplicate, default =1 ")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    if not inputbed:
        optparser.print_help()
        sys.exit(1)
    
    dupFeature(inputbed, options.outputfile, options.sequence, options.gcrange, options.flank)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)




