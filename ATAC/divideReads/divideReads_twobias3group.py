#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re,random
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
    biasATAC={}
    biasDNase={}
    biasShare={}
    inf = open(biasfile)
    for line in inf:
        ll = line.split()
        if ll[3] == "ATAC":
            biasATAC[ll[0]]= 0
        elif ll[3] == "DNase":
            biasDNase[ll[0]]= 0
        elif ll[3] == "share":
            biasShare[ll[0]]= 0
        else:
            print ll
    inf.close()
    return biasATAC,biasDNase,biasShare

def divideReadsFromBias(infile,outname,Biasfile,sequence,flank):
    genome = twobitreader.TwoBitFile(sequence)
    BA,BD,BS = dividebias_list(Biasfile)
    inf = open(infile)
    outf = open(outname+"_all.bed",'w')
    outfATAC = open(outname+"_biasATAC.bed",'w')
    outfDNase = open(outname+"_biasDNase.bed",'w')
    outfShare = open(outname+"_biasShare.bed",'w')
    random.seed(1228)

    for line in inf:
        ll = line.strip().split("\t")
        if ll[5] == '+' :
            seq = genome[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
        elif ll[5] == '-' :
            seq = rev(genome[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper())
        else:
            rdnum = random.randint(0,1)
            if rdnum == 0:
                seq = genome[ll[0]][(int(ll[1])-flank):(int(ll[1])+flank)].upper()
                newll = [ll[0],ll[1],str(int(ll[1])+1),ll[3],ll[4],"+"]
            else:
                seq = rev(genome[ll[0]][(int(ll[2])-flank):(int(ll[2])+flank)].upper())
                newll = [ll[0],str(int(ll[2])-1),ll[2],ll[3],ll[4],"-"]
            line= "\t".join(newll)+"\n"
                    
        if BA.has_key(seq):
            outfATAC.write(line)
        elif BD.has_key(seq):
            outfDNase.write(line)
        elif BS.has_key(seq):
            outfShare.write(line)
        
        if "N" in seq or len(seq) < 8:
            pass
        else:
            outf.write(line)
        
    outf.close()
    outfATAC.close()
    outfDNase.close()
    outfShare.close()
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
                         help="new format, seqtype,atacbias,dnasebias,group")              
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



