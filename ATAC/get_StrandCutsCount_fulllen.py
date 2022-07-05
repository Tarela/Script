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

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
def bwsig_pattern(bwfile,chrm,start,end,points):
    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
    sigPat = sp(cmd)[0].strip().split("\t")
    return sigPat


def get_signal(inputfile,output,plusbw,minusbw,readscut,bwfolder):
    if not bwfolder:
        bwfolder = ""
    if not bwfolder.endswith('/'):
        bwfolder += '/'
    #bwHs = []
    if not "/" in plusbw:
        plusbw = bwfolder + plusbw
    if not "/" in minusbw:
        minusbw = bwfolder + minusbw
    
    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        if "_" in ll[0]:
            continue
        #center = (int(ll[1]) + int(ll[2]))/2
        S = int(ll[1])# max(0,center - extend)
        E = int(ll[2])#center + extend
        result_p = bwsig_pattern(plusbw,ll[0],S,E,E-S)
        result_m = bwsig_pattern(minusbw,ll[0],S,E,E-S)
        result = result_p + result_m
        if 'n/a' in result or len(result) != (E-S)*2:
            continue#print result
        numpy_result = numpy.array(map(float,result))
        numpy_result[numpy_result > readscut] = readscut
        reads_count = sum(numpy_result)
        newll = ll + [reads_count]
        outf.write("\t".join(map(str,newll))+"\n")

    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """calculate the total reads count within 200bp of peak summits, reads inputed in plus_bw (1bp) and minus_bw (1bp)"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")
    optparser.add_option("-p","--plusbw",dest="plusbw",type="str",default = "",
                         help="")
    optparser.add_option("-n","--minusbw",dest="minusbw",type="str",default = "",
                         help="")
    optparser.add_option("--readscut",dest="readscut",type="int",default = "20",
                         help="if certain loci has more than 20 cuts on a strand, reduce it to 20")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",default = "",
                         help="")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    plusbw = options.plusbw
    minusbw = options.minusbw
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,plusbw,minusbw,options.readscut,options.bwfolder)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


