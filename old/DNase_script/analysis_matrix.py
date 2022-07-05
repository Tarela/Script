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

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def analysis(inputfile,out):
    outdict = {}
    for i in range(30):
        outdict[i] = [0]*100
    inf = open(inputfile)
    for line in inf:
        ll = line.split()
        obs = int(float(ll[0]))
        predict = int(float(ll[3]))
        if outdict.has_key(predict):
            if obs < 100:
                outdict[predict][obs] += 1
    inf.close()
    outf =open(out,'w')
    for i in sorted(outdict.keys()):
        newll = [i] + outdict[i]
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
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    out = options.out
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    analysis(inputfile,out)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


