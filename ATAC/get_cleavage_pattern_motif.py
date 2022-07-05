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
import string,numpy
from bx.bbi.bigwig_file import BigWigFile

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

#import subprocess
#def sp(cmd):
#    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
#    ac = a.communicate()
#    return ac
#def bwsig_pattern(bwfile,chrm,start,end,points):
#    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
#    sigPat = sp(cmd)[0].strip().split("\t")
#    return sigPat

def get_signal(inputfile,output,plus,minus,fulllen):

    plusbw = BigWigFile(open(plus, 'rb'))
    minusbw = BigWigFile(open(minus, 'rb'))

    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        motiflen = int(ll[2]) - int(ll[1])
        upstream_ext = fulllen/2 - motiflen/2
        try:
            if ll[5] == "+":
                start = int(ll[1])-upstream_ext
                end = start + fulllen
                forward_signal = list(plusbw.summarize(ll[0],start,end,end-start).sum_data)
                reverse_signal = list(minusbw.summarize(ll[0],start,end,end-start).sum_data)
            else:
                end = int(ll[2]) + upstream_ext
                start = end - fulllen
                forward_signal = list(minusbw.summarize(ll[0],start,end,end-start).sum_data)[::-1]
                reverse_signal = list(plusbw.summarize(ll[0],start,end,end-start).sum_data)[::-1]
        except:
            print ll
            forward_signal = [0]*(end-start)
            reverse_signal = [0]*(end-start)

        newll = ll + forward_signal + reverse_signal
        outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """calculate the average signal within full len of input peak, multiple bw files separated by comma"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="input bed file")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="output bedEX file")
    optparser.add_option("-p","--plusbw",dest="plusbw",type="str",default = "",
                         help="")
    optparser.add_option("-m","--minusbw",dest="minusbw",type="str",default = "",
                         help="")
    optparser.add_option("--fulllen",dest="fulllen",type="int",default = 100,
                         help="extend size from center")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    if not inputfile:
        optparser.print_help()
        sys.exit(1)

    get_signal(inputfile,output,options.plusbw,options.minusbw,options.fulllen)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

