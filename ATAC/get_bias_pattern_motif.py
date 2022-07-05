#!/usr/bin/env python
'''
Created on XXXX-XX-XX

@author: Tarela
'''
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging,time
import string,numpy

import twobitreader
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

#import subprocess
#def sp(cmd):
#    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
#    ac = a.communicate()
#    return ac
#def bwsig_pattern(bwfile,chrm,start,end,points):
#    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
#    sigPat = sp(cmd)[0].strip().split("\t")
#    return sigPat
def revcomp(seq):
    useseq = seq.upper()
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A' ,'N':'N'}
    seqt = list(useseq)
    seqt.reverse()
    r = ''.join( [ rc[x] for x in seqt] )
    return r

def readBG(bgmatrix):
    #Mgenomepos = 0
    #Nreadscount = 0
    pBG = {}
    inf = open(bgmatrix)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        name = ll[0]
#        pBG[name] = pow(numpy.e,float(ll[1]))
        pBG[name] = round(float(ll[2]),4)
    #    Mgenomepos += int(ll[3])
    #    Nreadscount += int(ll[2])        
    inf.close()
    return pBG,len(name)/2#, Nreadscount, Mgenomepos

def get_signal(inputfile,output,genome2bit,fulllen,biasmat):
    genome = twobitreader.TwoBitFile(genome2bit)
    biasdict,flank = readBG(biasmat)
    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        motiflen = int(ll[2]) - int(ll[1])
        upstream_ext = fulllen/2 - motiflen/2
        if float(ll[4]) <= 10000:
            continue
        if ll[5] == "+":
            start = int(ll[1])-upstream_ext
            end = start + fulllen

#            forward_signal = list(plusbw.summarize(ll[0],start,end,end-start).sum_data)
#            reverse_signal = list(minusbw.summarize(ll[0],start,end,end-start).sum_data)
        else:
            end = int(ll[2]) + upstream_ext
            start = end - fulllen
#            forward_signal = list(minusbw.summarize(ll[0],start,end,end-start).sum_data)[::-1]
#            reverse_signal = list(plusbw.summarize(ll[0],start,end,end-start).sum_data)[::-1]

        plusSequence = genome[chrm][(start-flank):(end+flank)].upper()
        minusSequence = genome[chrm][(start-flank+1):(end+flank+1)].upper()
        
        Plus_tmp = []
        Minus_tmp = []

        for i in range(fulllen):
            pseq = plusSequence[i:(i+2*flank)].upper()
            mseq = revcomp(minusSequence[i:(i+2*flank)]).upper()
            if biasdict.has_key(pseq):
                Plus_tmp.append(biasdict[pseq])
            else:
                Plus_tmp.append(0)
            if biasdict.has_key(mseq):
                Minus_tmp.append(biasdict[mseq])
            else:
                Minus_tmp.append(0)

        if ll[5] == "+":
            forward_signal = Plus_tmp
            reverse_signal = Minus_tmp
        else:
            forward_signal = Minus_tmp[::-1]
            reverse_signal = Plus_tmp[::-1]

        newll = ll + forward_signal + reverse_signal
        outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """calculate the average signal within full len of input peak, multiple bw files separated by comma"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="input bed file")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="output bedEX file")
    optparser.add_option("--fulllen",dest="fulllen",type="int",default = 100,
                         help="extend size from center")
    optparser.add_option("-s","--species",dest="species",type="str",default = "/scratch/sh8tv/Data/Genome/hg38/hg38.2bit",
                         help="genome 2bit sequence file")
    optparser.add_option("-b","--biasmat",dest="biasmat",type="str",
                         help="cutting bias matrix")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    if not inputfile:
        optparser.print_help()
        sys.exit(1)

    get_signal(inputfile,output,options.species,options.fulllen,options.biasmat)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

