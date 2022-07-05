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
import numpy
#from numpy import linalg as la

#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def sitepro_scan(peak,outname,w_plus,w_minus,Cspan):

    
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    outf_propPlus = open(outname + "_propcutPlus.bdg",'w')
    outf_propMinus = open(outname + "_propcutMinus.bdg",'w')

    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        if start-Cspan < 0:
            print ll
            continue
        plus_obj = w_plus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        minus_obj = w_minus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        if not plus_obj  :
            plus_vector = numpy.array([0]*(end-start+2*Cspan))+1
        else:
            plus_vector = plus_obj.sum_data + 1
        if not minus_obj  :
            minus_vector = numpy.array([0]*(end-start+2*Cspan))+1
        else:
            minus_vector = minus_obj.sum_data + 1
   
        roundN = 4 
        #### assign bias to bp and proportion
        for outpos in range(Cspan,(end-start+Cspan)):
                        
            this_plus_cuts_prop = round(plus_vector[outpos]/sum(plus_vector[(outpos-Cspan):(outpos+Cspan)]),roundN)
            this_minus_cuts_prop = round(minus_vector[outpos]/sum(minus_vector[(outpos-Cspan):(outpos+Cspan)]),roundN)
                    
            out_chrm = chrm
            out_start = start + outpos - Cspan 
            out_end = out_start+1
            
            outf_propPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cuts_prop]))+"\n")
            outf_propMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cuts_prop]))+"\n")

    outf_propPlus.close()
    outf_propMinus.close()

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
                         help="5 column summit200bp file ,sorted")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="name of output file")
    optparser.add_option("-p","--plusbw",dest="plusbw",type="str",
                         help="plus cuts in bw format")
    optparser.add_option("-n","--minusbw",dest="minusbw",type="str",
                         help="minus cuts in bw format")

#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 25,
                         help="region for get total signal in single bp, default = 25 means +-25bp(total 50bp) signal as total for each bp")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    outname = options.outname
    plusbw = options.plusbw
    minusbw = options.minusbw
    Cspan = options.Cspan
    sitepro_scan(interval,outname,plusbw,minusbw,Cspan)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

