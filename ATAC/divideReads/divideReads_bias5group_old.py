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

def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
  #      print seq
        nmer_seq[seq] = 0
    del allseq
    return nmer_seq
    
def dividebias_list(biasfile):
    biaslist = []
    bias1=[]
    bias2=[]
    bias3=[]
    bias4=[]
    bias5=[]
    inf = open(biasfile)
    for line in inf:
        ll = line.split()
        thisll = [ll[0],float(ll[2])]
        biaslist.append(thisll)
    inf.close()
    sortbias = sorted(biaslist,key=lambda x:x[1])
    Neachgroup = len(sortbias)/5
    for idx in range(len(sortbias)):
        if idx < Neachgroup:
            bias1.append(sortbias[idx][0])
        elif idx < Neachgroup*2:
            bias2.append(sortbias[idx][0])
        elif idx < Neachgroup*3:
            bias3.append(sortbias[idx][0])
        elif idx < Neachgroup*4:
            bias4.append(sortbias[idx][0])
        else:
            bias5.append(sortbias[idx][0])
    return bias1,bias2,bias3,bias4,bias5

def divideReadsFromBias(infile,outname,Biasfile,sequence,flank):
    genome = twobitreader.TwoBitFile(sequence)
    B1,B2,B3,B4,B5 = dividebias_list(Biasfile)
    inf = open(infile)
    outf1 = open(outname+"_biasG1.bed",'w')
    outf2 = open(outname+"_biasG2.bed",'w')
    outf3 = open(outname+"_biasG3.bed",'w')
    outf4 = open(outname+"_biasG4.bed",'w')
    outf5 = open(outname+"_biasG5.bed",'w')
    
    for line in inf:
        ll = line.strip().split("\t")
        if ll[5] == '+' :
            seq = genome[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
        elif ll[5] == '-' :
            seq = rev(genome[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper())
        else:
            print ll
                    
        if seq in B1:
            outf1.write(line)
        elif seq in B2:
            outf2.write(line)
        elif seq in B3:
            outf3.write(line)
        elif seq in B4:
            outf4.write(line)
        elif seq in B5:
            outf5.write(line)
        elif "N" in seq or len(seq) < 8:
            pass
        else:
            print seq,ll

    outf1.close()
    outf2.close()
    outf3.close()
    outf4.close()
    outf5.close()
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
    optparser.add_option("-i","--inputfile",dest="infile",type="str",
                         help="ATAC/DNase align bed file")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="")              
    optparser.add_option("-b","--biasfile",dest="biasfile",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    infile = options.infile
    biasfile = options.biasfile
    out = options.outname
    seq = options.sequence
    flank = options.flank
    if not infile or not biasfile:
        optparser.print_help()
        sys.exit(1)
    
    divideReadsFromBias(infile,out,biasfile,seq,flank)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



