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

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def duplication(inputfile,outfile):
    inf = open(inputfile)
    last = inf.readline()
    outf = open(outfile,'w')
    for line in inf:
        if line == last:
            pass
        else:
            outf.write(line)
        last = line
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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--outfile",dest="outfile",type="str",default = "",
                         help="")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outfile = options.outfile
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    duplication(inputfile,outfile)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

