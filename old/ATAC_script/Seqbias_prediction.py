'''
Created on XXXX-XX-XX

@author: Tarela
'''
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
import math,time
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader

#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def readBG(bgmatrix):
    pBG = {}
    nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        ll = line.split()
        name = ll[0]
        pBG[name] = float(ll[1])
        nBG[name] = float(ll[2])
    inf.close()
    return pBG,nBG
def sitepro_scan(peak,out_bed,w_plus,w_minus,bgmatrix,span,gen,lflank,rflank,offset,bpshift):
    nmer = lflank + rflank
    genome = twobitreader.TwoBitFile(gen)
    pBG,nBG = readBG(bgmatrix)
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    outf = open(out_bed,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        ## remove overflow
        if start - span -lflank - offset <= 0:
            print line
            continue
        ## get cleavage
        p_sum = list(w_plus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        n_sum = list(w_minus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        ## get seqbias

        if 'N' in genome[chrm][min((start-span+1-offset +bpshift-lflank),(start-span+1 -bpshift-rflank) ):max((end+span+offset+lflank-bpshift),(end+span + bpshift + rflank))].upper():
            print line
            print genome[chrm][min((start-span+1-offset +bpshift-lflank),(start-span+1 -bpshift-rflank) ):max((end+span+offset+lflank-bpshift),(end+span + bpshift + rflank))].upper()
            continue
        praw=[]
        nraw=[]
        px = []
        nx = []
        psqx=[]
        nsqx=[]
        pmax = []
        nmax = []
        
        for bp1 in range(-span,end-start+span):
            loci = start + bp1
            pseq = genome[chrm][(loci + bpshift - lflank) : (loci + bpshift + rflank)].upper()
            pseq_apart = genome[chrm][(loci+offset -bpshift-rflank):(loci+offset -bpshift+lflank)].upper()
            nseq = genome[chrm][(loci+1 -bpshift-rflank) : (loci+1 -bpshift+lflank)].upper()
            nseq_apart = genome[chrm][(loci+1-offset +bpshift-lflank):(loci+1-offset +bpshift+rflank)].upper()
            praw.append(pBG[pseq])
            nraw.append(nBG[nseq])
            px.append(pBG[pseq] * nBG[pseq_apart])
            
            nx.append(nBG[nseq] * pBG[nseq_apart])
            psqx.append(math.sqrt(pBG[pseq] * nBG[pseq_apart]))
            nsqx.append(math.sqrt(nBG[nseq] * pBG[nseq_apart]))
            pmax.append(  max(pBG[pseq] , nBG[pseq_apart]) )
            nmax.append(  max(nBG[nseq] , pBG[nseq_apart]) )
            
        ## get predicted seqbias
        praw_assign =[]
        nraw_assign =[]
        px_assign  = []
        nx_assign  = []
        psqx_assign =[]
        nsqx_assign =[]
        pmax_assign  = []
        nmax_assign  = []

        for bp in range(len(p_sum)- 2*span):
            
            ptotal = sum(p_sum[bp:(bp+2*span)])*1.0
            ntotal = sum(n_sum[bp:(bp+2*span)])*1.0
            
            praw_assign.append(ptotal * praw[bp+span]/sum(praw[bp:(bp+2*span)]))
            nraw_assign.append(ntotal * nraw[bp+span]/sum(nraw[bp:(bp+2*span)]))
            px_assign.append(ptotal * px[bp+span]/sum(px[bp:(bp+2*span)]))
            nx_assign.append(ntotal * nx[bp+span]/sum(nx[bp:(bp+2*span)]))
            psqx_assign.append(ptotal * psqx[bp+span]/sum(psqx[bp:(bp+2*span)]))
            nsqx_assign.append(ntotal * nsqx[bp+span]/sum(nsqx[bp:(bp+2*span)]))
            pmax_assign.append(ptotal * pmax[bp+span]/sum(pmax[bp:(bp+2*span)]))
            nmax_assign.append(ntotal * nmax[bp+span]/sum(nmax[bp:(bp+2*span)]))
 


        ### write  real cleavage , seqbias , seqbias predicted cleavage
        newll = ll + p_sum[span:(len(p_sum)-span)] + n_sum[span:(len(n_sum)-span)] + praw[span:(len(praw)-span)] + nraw[span:(len(nraw)-span)] + praw_assign + nraw_assign + px_assign + nx_assign + psqx_assign + nsqx_assign + pmax_assign + nmax_assign
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
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="if raw , -i is 5 column summit200bp file")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/Bam/GM12878_ATACseq_50k_p.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/Bam/GM12878_ATACseq_50k_m.bw",
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/biasoffset/up3down3/GM12878_ATACseq_50k_bias_onPeak_shift0.txt",
                         help="sequence bias matrix")
#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 25,
                         help="region for get total signal in single bp, default = 25 means +-25bp(total 50bp) signal as total for each bp")

    optparser.add_option("--genome",dest="genome",type="str",default = "/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit",
                         help="2bit format")
    optparser.add_option("--left",dest="leftflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--right",dest="rightflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--offset",dest="offset",type="int",default = 9,
                         help="offset related, distance of pair of +/- related cut,default = 9")
    optparser.add_option("--bpshift",dest="bpshift",type="int",default = 0,
                         help="bp of shift center when calculate bias for given cut site,default is 0 ")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    out_bed = options.output
    w_plus = options.w_plus
    w_minus = options.w_minus
    bgmatrix = options.bgmatrix
    gen = options.genome
    lflank = options.leftflank
    rflank = options.rightflank
    Cspan = options.Cspan
    offset = options.offset
    sitepro_scan(interval,out_bed,w_plus,w_minus,bgmatrix,Cspan,gen,lflank,rflank,offset,options.bpshift)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

