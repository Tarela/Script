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
import twobitreader,numpy,random
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
    
def readDict(BFile):
    outD = {}
    outlist = []
    inf = open(BFile)
    for line in inf:
        ll =line.split()
        outD[ll[0]] = float(ll[2])
        outlist.append(float(ll[2]))
    inf.close()
    return outD,numpy.median(outlist)

def divideReadsFromBias(infile,outname,BiasfileA,BiasfileD,sequence,flank):
    genome = twobitreader.TwoBitFile(sequence)
    biasATAC,medATAC = readDict(BiasfileA)
    biasDNase,medDNase = readDict(BiasfileD)
    
    inf = open(infile)
    outf = open(outname+"_all.bed",'w')
    outf11 = open(outname+"_biasG11.bed",'w')
    outf10 = open(outname+"_biasG10.bed",'w')
    outf01 = open(outname+"_biasG01.bed",'w')
    outf00 = open(outname+"_biasG00.bed",'w')
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
        
        if  not biasATAC.has_key(seq) or not biasDNase.has_key(seq):
            continue
            
        biasA = biasATAC[seq]
        biasD = biasDNase[seq]
        outf.write(line)
        if biasA > medATAC:
            if biasD > medDNase:
                outf11.write(line)
            else:
                outf10.write(line)
        else:
            if biasD > medDNase:
                outf01.write(line)
            else:
                outf00.write(line)

    outf11.close()
    outf10.close()
    outf01.close()
    outf00.close()
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
    optparser.add_option("-i","--inputfile",dest="infile",type="str",
                         help="ATAC/DNase align bed file")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="")              
    optparser.add_option("-a","--biasfileATAC",dest="biasfileATAC",type="str",
                         help="")              
    optparser.add_option("-d","--biasfileDNase",dest="biasfileDNase",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    infile = options.infile
    out = options.outname
    seq = options.sequence
    flank = options.flank
    if not infile or not out:
        optparser.print_help()
        sys.exit(1)
    
    divideReadsFromBias(infile,out,options.biasfileATAC,options.biasfileDNase,seq,flank)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



