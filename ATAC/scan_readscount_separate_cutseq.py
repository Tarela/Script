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
import math,time,random
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader
#import scipy.stats.distributions
from copy import deepcopy
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def revcomp(seq):
    useseq = seq.upper()
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A' ,'N':'N'}
    seqt = list(useseq)
    seqt.reverse()
    r = ''.join( [ rc[x] for x in seqt] )
    return r


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

def make_rd_dict(n):
    Nkeys = 4**n
    nmer_seq = {}
    for i in range(1,(Nkeys+1)):
        nmer_seq['rd'+str(i)]=0
    return nmer_seq

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

def get_regionLevel_reads(inbed,outputname,plusbw,minusbw,species,flank,ext):
    genome = twobitreader.TwoBitFile("/scratch/sh8tv/Data/Genome/%s/%s.2bit"%(species,species))
    countdict_template = make_nmer_dict(2*flank)
    rddict_template = make_rd_dict(2*flank)

    plusBWH = BigWigFile(open(plusbw, 'rb'))
    minusBWH = BigWigFile(open(minusbw, 'rb'))

    random.seed(1228)

    inf = open(inbed)
    outf = open(outputname + "_seqtype.bed",'w')
    outfRD = open(outputname + "_rd.bed",'w')
    seqtypes = sorted(countdict_template.keys())
    newll_seq = []
    newll_rd = []
    infileLen = len(inf.readline().split())
    for i in range(infileLen):
        newll_seq.append("C"+str(i))
        newll_rd.append("C"+str(i))
    newll_seq += sorted(countdict_template.keys())
    newll_rd += sorted(rddict_template.keys())

    outf.write("\t".join(newll_seq)+"\n")
    outfRD.write("\t".join(newll_rd)+"\n")

    inf.seek(0)
    for line in inf:
        Sdict = deepcopy(countdict_template)
        Rdict = deepcopy(rddict_template)
        ll = line.split()
        chrm = ll[0]
        if "_" in chrm:
            continue
        center = (int(ll[1]) + int(ll[2]))/2
        start = max(0,center-ext)
        end = center + ext
        plusSig_raw = plusBWH.summarize(ll[0],start,end,end-start)#.sum_data
        minusSig_raw = minusBWH.summarize(ll[0],start,end,end-start)#.sum_data
        if type(plusSig_raw) == None or type(minusSig_raw) == None:
            continue
        plusSig = plusSig_raw.sum_data
        minusSig = minusSig_raw.sum_data
        plusSequence = genome[chrm][(start-flank):(end+flank)].upper()
        minusSequence = genome[chrm][(start-flank+1):(end+flank+1)].upper()

        for i in range(len(plusSig)):
            #position = start + i
            pcuts = plusSig[i]
            if pcuts > 0:
                pseq = plusSequence[i:(i+2*flank)].upper()
                #pseqRV = revcomp(plusSequence_reverse[i:(i+2*flank)]).upper()
                if not "N" in pseq :#and not 'N' in pseqRV:
                #    p_out = seq2biasParm(pseq,B,simplex_code)
                #    plus_data += pcuts*p_out
                    Sdict[pseq] += 1#pcuts
                    Rdict["rd"+str(random.randint(1,4**(2*flank)))] += 1#pcuts
                    #plus_readscount += pcuts
                    #plus_biassum += biasdict[pseq]*pcuts
                    #plus_biasCB += (biasdict[pseq]+biasdict[pseqRV] ) *pcuts/2

                    #print i,pcuts,plus_readscount   
        for i in range(len(minusSig)):
            #position = start + i
            mcuts = minusSig[i]
            if mcuts > 0:
#                tmpseq = minusSequence[i:(i+2*flank)]
                mseq = revcomp(minusSequence[i:(i+2*flank)]).upper()
                #mseqRV = minusSequence_reverse[i:(i+2*flank)].upper()
                if not "N" in mseq :#and not "N" in mseqRV:
                #    m_out = seq2biasParm(mseq,B,simplex_code)
                #    minus_data += mcuts*m_out
                    Sdict[mseq] += 1#mcuts
                    Rdict["rd"+str(random.randint(1,4**(2*flank)))] += 1#mcuts
                    #minus_readscount += mcuts
                    #minus_biassum += biasdict[mseq]*mcuts
                    #minus_biasCB += (biasdict[mseq]+biasdict[mseqRV] ) *mcuts/2
#                print chrm,start,end,i,mcuts,minus_biassum,minus_biasCB
        #plus_biasave = plus_biassum / plus_readscount
        #minus_biasave = minus_biassum / minus_readscount
        #newll = ll + [plus_readscount,minus_readscount,plus_biassum,minus_biassum]#plus_biassum,minus_biassum,plus_biasCB,minus_biasCB] #+ list(plus_data) + list(minus_data)
        newll_seq = ll + [Sdict[x] for x in sorted(Sdict.keys())]
        newll_rd = ll + [Rdict[x] for x in sorted(Rdict.keys())]
        outf.write("\t".join(map(str,newll_seq))+"\n")
        outfRD.write("\t".join(map(str,newll_rd))+"\n")

    inf.close()
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
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="inputbed file")
    optparser.add_option("-o","--outputname",dest="outputname",type="str",
                         help="output bed file, reads count for different cleavage sequence")              
    optparser.add_option("-e","--ext",dest="ext",type="int",default=1000,
                         help="extend size from the center of input bed region")              
    optparser.add_option("-s","--species",dest="species",type="str",default='hg38',
                         help="genome version, choose from hg38, mm10, sacCer3")
    optparser.add_option("-p","--plusbw",dest="plusbw",type="str",
                         help="1bp resolution cleavage bigwig for plus strand cuts")
    optparser.add_option("-m","--minusbw",dest="minusbw",type="str",
                         help="1bp resolution cleavage bigwig for minus strand cuts")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputname = options.outputname
    plusbw = options.plusbw
    minusbw = options.minusbw

    if not inputbed or not outputname:
        optparser.print_help()
        sys.exit(1)
    get_regionLevel_reads(inputbed,outputname,plusbw,minusbw,options.species,options.flank,options.ext)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)




