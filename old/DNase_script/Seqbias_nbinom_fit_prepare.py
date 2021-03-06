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

import scipy.stats.distributions

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
def sitepro_scan(peak,outp,outn,w_plus,w_minus,bgmatrix,span,gen,lflank,rflank):
    nmer = lflank + rflank
    genome = twobitreader.TwoBitFile(gen)
    pBG,nBG = readBG(bgmatrix)
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    outfp = open(outp,'w')
    outfn = open(outn,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        ## remove overflow
        if start - span -lflank <= 0:
            continue
        ## get cleavage
        p_sum = list(w_plus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        n_sum = list(w_minus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        ## get seqbias
        seq = genome[chrm][(start - span - lflank):(end + span + rflank)]
        if 'N' in seq.upper():
            continue
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]### bias
        for k in range(len(pseq)  +1 - nmer):
            p.append(pBG[pseq[k:(k+nmer)].upper()])
            n.append(nBG[nseq[k:(k+nmer)].upper()])

        for bp in range(len(p_sum)- 2*span):
            ptotal = sum(p_sum[bp:(bp+2*span)])### total
            ntotal = sum(n_sum[bp:(bp+2*span)])
            pc = int(p_sum[bp+span])#### observation cut
            nc = int(n_sum[bp+span])            
            pbias = p[bp+span]
            nbias = n[bp+span]
            pbgtotal = sum(p[bp:(bp+span*2)])
            nbgtotal = sum(n[bp:(bp+span*2)])
            paraw = (pbias/pbgtotal)*ptotal
            naraw = (nbias/nbgtotal)*ntotal

            outfp.write( "\t".join(map(str,[pc,ptotal,pbias,paraw])) + "\n" )
            outfn.write( "\t".join(map(str,[nc,ntotal,nbias,paraw])) + "\n" )
    outfp.close()
    outfn.close()
    inf.close()
            ## calculate real cut distribution



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
                         help="-i is 5 column summit200bp file")
    optparser.add_option("-n","--outname",dest="name",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_plus_50_100.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_minus_50_100.bw",
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="sequence bias matrix")
#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 25,
                         help="region for get total signal in single bp, default = 25 means +-25bp(total 50bp) signal as total for each bp")
    optparser.add_option("--genome",dest="genome",type="str",
                         help="2bit format")
    optparser.add_option("--left",dest="leftflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--right",dest="rightflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)


    interval = options.interval
    outp = options.name + '_plus_matrix.txt'
    outn = options.name + '_minus_matrix.txt'
    w_plus = options.w_plus
    w_minus = options.w_minus
    bgmatrix = options.bgmatrix
    gen = options.genome
    lflank = options.leftflank
    rflank = options.rightflank
    
    sitepro_scan(interval,outp,outn,w_plus,w_minus,bgmatrix,options.Cspan,gen,lflank,rflank)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

