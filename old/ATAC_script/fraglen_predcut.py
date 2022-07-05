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
import numpy
#from cistrome import regions as c
import twobitreader
from bx.bbi.bigwig_file import BigWigFile
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def fragment_v_predcut(pefrag,output,bwp,bwn):

    inf = open(pefrag)
    outf = open(output,'w')
    pH = BigWigFile(open(bwp, 'rb'))
    mH = BigWigFile(open(bwn, 'rb'))
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        pcut = int(ll[1])
        ncut = int(ll[2])-1
        fraglen = ncut-pcut+1
        pPred = float(pH.summarize(chrm,pcut,pcut+1,1).sum_data)
        nPred = float(mH.summarize(chrm,ncut,ncut+1,1).sum_data)
        if pPred == -1 or nPred == -1:
            continue
        newll = [chrm,pcut,ncut,fraglen,pPred,nPred]
        outf.write("\t".join(map(str,newll))+"\n")
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
    optparser.add_option("-i","--pefrag",dest="pefrag",type="str",
                         help="PEfragment bed")
    optparser.add_option("-o","--output",dest="out",type="str",
                         help="")
    optparser.add_option("-p","--predbwplus",dest="wp",type="str",default = '/mnt/Storage/home/huse/Project/ATAC/Result/simplex_encode/FPprediction/high_bias_cause/Data/ATAC_u4d4s0_sorted_encoding_9bpbias_plus.bw',
                         help="seqbias pred cut bigWiggle file")
    optparser.add_option("-n","--predbwminus",dest="wm",type="str",default = '/mnt/Storage/home/huse/Project/ATAC/Result/simplex_encode/FPprediction/high_bias_cause/Data/ATAC_u4d4s0_sorted_encoding_9bpbias_minus.bw',
                         help="seqbias pred cut bigWiggle file")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    pefrag = options.pefrag
    out = options.out
    w_plus = options.wp
    w_minus = options.wm
    if not pefrag:
        optparser.print_help()
        sys.exit(1)
    
    fragment_v_predcut(pefrag,out,w_plus,w_minus)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)




