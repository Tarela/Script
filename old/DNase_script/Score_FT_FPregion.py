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
import math
#from CistromeAP.taolib.CoreLib.BasicStat.Func import *
#from CistromeAP.jianlib.BwReader import BwIO
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

    
def getsignal(inputfile,outputfile,pcut,pspan,FPregion):

    
#    p=BwIO(pcut)
#    chrom_len = {}
#    for i in p.chromosomeTree['nodes']:
#        chrom_len[i['key']] = i['chromSize']
    pcutbw = BigWigFile(open(pcut, 'rb'))
    FPbw = BigWigFile(open(FPregion,'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    inf.seek(0)
    outf = open(outputfile,'w')

    for line in inf:
        ll = line.split()
#        if not chrom_len.has_key(ll[0]):
#            continue
        cut = list(pcutbw.summarize(ll[0],int(ll[1]) + ml/2 -pspan ,int(ll[1]) + ml/2 +pspan ,2*pspan).sum_data)
        TC = sum(cut)
        C = sum(cut[(pspan-ml/2) : (pspan-ml/2+ml)])
        L = sum(cut[(pspan-ml/2-ml):(pspan-ml/2)])
        R = sum(cut[(pspan-ml/2+ml):(pspan-ml/2+2*ml)])
        FOS = -1*( (C+1)/(R+1) + (C+1)/(L+1) )
        try:
            FP_bw = map(float,list(FPbw.summarize(ll[0],int(ll[1]) ,int(ll[2]),int(ll[2])-int(ll[1])).sum_data))
        except:
            FP_bw = [0.0]*(int(ll[2])-int(ll[1]))
        minFPbw = min(FP_bw)
        maxFPbw = max(FP_bw)
        newll = ll + [TC,FOS,minFPbw,maxFPbw]
        outf.write("\t".join(map(str,newll))+"\n")

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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")                         

    optparser.add_option("-w","--cutbw",dest="cut",type="str",
                         help="")
    optparser.add_option("-f","--fpregion",dest="fpregion",type="str",
                         help="neph's FPregion (k562)")
    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 100,
                         help="window size for tagcount and dymDHS, default = 100 , means total length = 200")

 
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    cutbw = options.cut
    fpregion = options.fpregion
    pspan = options.pspan
    
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,cutbw,pspan,fpregion)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
