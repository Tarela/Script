#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description: this one seems to have bugs in fetch sequence

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
import numpy
from numpy import linalg as la
from copy import deepcopy
#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def peak_cut_count(peak,out_bed,w_plus,w_minus,span,gen):
 
    genome = twobitreader.TwoBitFile(gen)
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    outf = open(out_bed,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        p_sum = float(w_plus_H.summarize(chrm,start-span,end+span,1).sum_data)
        n_sum = float(w_minus_H.summarize(chrm,start-span,end+span,1).sum_data)
        newll = ll+[p_sum+n_sum]
        outf.write("\t".join(map(str,newll))+"\n")
    #print "predict cut time :",time.time()-t

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
    optparser.add_option("-s","--span",dest="span",type="int",default=200,
                         help="")
    optparser.add_option("--genome",dest="genome",type="str",default = "/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit",
                         help="2bit format")

    (options,args) = optparser.parse_args()

    if not options.interval or not options.output:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    out_bed = options.output
    w_plus = options.w_plus
    w_minus = options.w_minus
    peak_cut_count(interval,out_bed,w_plus,w_minus,options.span,options.genome)
    
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

