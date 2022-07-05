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
import string
from CistromeAP.taolib.CoreLib.BasicStat.Func import *
from CistromeAP.jianlib.BwReader import BwIO
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

def get_signal(inputfile,output,signalbw,extend):
    signalbw = signalbw.strip().strip(',').split(',')
    
    p=BwIO(signalbw[0])
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    bwHandle = []
    for k in signalbw:
        bwHandle.append(BigWigFile(open(k, 'rb')))
    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        inputlen = len(ll)
        if not chrom_len.has_key(ll[0]):
            continue
        for bwH in bwHandle:
            S = (int(ll[1]) + int(ll[2]))/2
            E = (int(ll[1]) + int(ll[2]))/2 + 1
            try:
                signal=bwH.summarize(ll[0],max(0,S-extend),E+extend,1)
            except:
                break
            if float(signal.valid_count) == 0:
                ll.append('0')
            else:
                ll.append(str(float(signal.sum_data/signal.valid_count)))
        if len(ll) == ( inputlen + len(bwHandle)  ):
            outf.write("\t".join(ll)+"\n")
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
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")
    optparser.add_option("-w","--signalbw",dest="bw",type="str",default = "",
                         help="")
    optparser.add_option("--ext",dest="ex",type="int",default = 0,
                         help="")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    signalbw = options.bw
    extend = options.ex
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,signalbw,extend)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

