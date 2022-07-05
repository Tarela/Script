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
def SOBfetchSEQ(fullseq,seqlen):
    if seqlen == 11:
        outseq = fullseq[:3]+fullseq[4]+fullseq[8:10]+fullseq[11:15]+fullseq[16]
    elif seqlen == 10:
        outseq = fullseq[:3]+fullseq[4]+fullseq[8:10]+fullseq[11:14]+fullseq[16]
    elif seqlen == 8:
        outseq = fullseq[0]+fullseq[2]+fullseq[8:10]+fullseq[11:14]+fullseq[16]
    elif seqlen == 6:
        outseq = fullseq[0]+fullseq[2]+fullseq[8:10]+fullseq[11:13]
    elif seqlen == 4:
        outseq = fullseq[2]+fullseq[8]+fullseq[11:13]
    else:
        outseq = "NA"
    return outseq

def seqbias(peak,tag,out,sequence,kmer):
    genome = twobitreader.TwoBitFile(sequence) 
    pcut = make_nmer_dict(kmer)
    bgseq = make_nmer_dict(kmer)
    inf = open(peak)
    for line in inf:
        ll = line.strip().split("\t")
        seq = genome[ll[0]][int(ll[1]):int(ll[2])]

        for i in range(len(seq)):
            subseq_template = seq[(i-6):(i+11)]
            if len(subseq_template) != 17:
                continue
            subseq = SOBfetchSEQ(subseq_template,kmer)
            if bgseq.has_key(subseq):
                bgseq[subseq]+=1
            else:
                pass
    inf.close()
    inf = open(tag)
    PEtag = 0
    for line in inf:
        ll = line.strip().split("\t")
        if len(ll) < 6:
            PEtag = 1
        if PEtag == 1 or ll[5] == '+' or ll[5] == "." :
            rawseq = genome[ll[0]][(int(ll[1])-6):(int(ll[1])+11)].upper()
            if len(rawseq) != 17:
                continue
            seq = SOBfetchSEQ(rawseq,kmer)
            if pcut.has_key(seq):
                pcut[seq]+=1
            else:
                #print seq
                pass
        if PEtag == 1 or ll[5] == '-' or ll[5] == ".":
            rawseq = rev(genome[ll[0]][(int(ll[2])-11):(int(ll[2])+6)].upper())
            if len(rawseq) != 17:
                continue
            seq = SOBfetchSEQ(rawseq,kmer)
            if pcut.has_key(seq):
                pcut[seq]+=1
            else:
                #print seq
                pass
        
    inf.close()
    outf = open(out,'w')
    for seqtype in sorted(pcut.keys()):
        if bgseq[seqtype] == 0:
            pbias = -1
        else:  
            pbias = float(pcut[seqtype])/float(bgseq[seqtype])
        #nbias = float(ncut[seqtype])/float(bgseq[seqtype])
        #outf.write("\t".join(map(str,[seqtype,pcut[seqtype]]))+"\n")
        outf.write("\t".join(map(str,[seqtype,pbias,pcut[seqtype],bgseq[seqtype]]))+"\n")
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-p","--peak",dest="peak",type="str",
                         help="peak region for calculate background (k-mer occurancy)")
    optparser.add_option("-t","--tag",dest="tag",type="str",
                         help="6-column reads file (bed) for calculate foreground (cuts occyrancy)")              
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-k","--kmer",dest="kmer",type="int",default=8,
                         help="kmer length , default =8")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    peak = options.peak
    tag = options.tag
    out = options.out
    seq = options.sequence
    kmer = options.kmer
    if not tag or not peak:
        optparser.print_help()
        sys.exit(1)
    if not kmer in [4,6,8,10,11]:
        print "kmer choose from 4,6,8,10,11"
        optparser.print_help()
        sys.exit(1)
    

    seqbias(peak,tag,out,seq,kmer)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



