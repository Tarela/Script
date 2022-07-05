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

#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def bp2number(inputfile,outputfile):
    inf =open(inputfile)
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()
        seq = ll[0]
        outnum =[]
        for bp in seq:
            if bp =="A":
                outnum.append(1)
            elif bp == "C":
                outnum.append(2)
            elif bp == "G":
                outnum.append(3)
            elif bp == "T":
                outnum.append(4)
        newll = ll + outnum
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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")


    (options,args) = optparser.parse_args()

    if not options.inputfile:
        optparser.print_help()
        sys.exit(1)

    inputfile = options.inputfile
    output = options.output
    
    bp2number(inputfile,output)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

