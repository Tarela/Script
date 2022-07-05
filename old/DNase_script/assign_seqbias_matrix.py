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
def sitepro_scan(peak,out,w_plus,w_minus,bgmatrix,span,gen,lflank,rflank):
    nmer = lflank + rflank
    genome = twobitreader.TwoBitFile(gen)
    pBG,nBG = readBG(bgmatrix)
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf = open(out,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        if start - span -lflank <= 0:
            continue
        p_sum = list(w_plus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        n_sum = list(w_minus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        seq = genome[chrm][(start - span - lflank):(end + span + rflank)]
        if 'N' in seq.upper():
            continue
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - nmer):
            p.append(pBG[pseq[k:(k+nmer)].upper()])
            n.append(nBG[nseq[k:(k+nmer)].upper()])
        p_assign = []
        n_assign = []
        for bp in range(len(p_sum)- 2*span):
            ptotal = sum(p_sum[bp:(bp+2*span)])
            ntotal = sum(n_sum[bp:(bp+2*span)])
            pbias_per = p[bp+span]*1.0/sum(p[bp:(bp+2*span)])
            nbias_per = n[bp+span]*1.0/sum(n[bp:(bp+2*span)])
            p_assign.append(pbias_per*ptotal)
            n_assign.append(nbias_per*ntotal)
            newll = [p_sum[(bp+span)],ptotal,p[bp+span],p_assign]
            outf.write("\t".join(map(str,newll))+"\n")
            newll = [n_sum[(bp+span)],ntotal,n[bp+span],n_assign]
            outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
        #print len(p_assign)    
            
        #print len(p_sum[span:(len(p_sum)-span)]),len(n_sum[span:(len(n_sum)-span)])
        #print len(seq)  
        #print len(p[span:(len(p)-span)])     
        #sys.exit(1)
        #fp=(ll+p_sum+m_sum)
        #newline = "\t".join(map(str,fp))+"\n"
        #outf.write(newline)
    
    #outf.close()
              #  print score

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
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_plus_50_100.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_minus_50_100.bw",
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="sequence bias matrix")
                                                                           
#========minor options=============
    optparser.add_option("-s","--span",dest="span",type="int",default = 10,
                         help="region for get total signal in single bp, default = 10 means +-10bp(total 20bp) signal as total for each bp")
    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
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
    output = options.output
    w_plus = options.w_plus
    w_minus = options.w_minus
    bgmatrix = options.bgmatrix
    span = options.span
    gen = options.genome
    lflank = options.leftflank
    rflank = options.rightflank
    
    sitepro_scan(interval,output,w_plus,w_minus,bgmatrix,span,gen,lflank,rflank)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

