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

def fraglen(inputfile,output):
    inf = open(inputfile)
    D = {}
    for i in range(600):
        D[i]=0
    for line in inf:
        if D.has_key(int(line.split()[3])):
            D[int(line.split()[3])]+=1    
    inf.close()
    length = []
    times = []
    for i in sorted(D.keys()):
        length.append(i)
        times.append(D[i])
    outf = open(output+'.r','w')
    outf.write('fraglen<-c('+str(length)[1:-1]+')\n')
    outf.write('%s<-c('%(output)+str(times)[1:-1]+')\n')
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
    optparser.add_option("-o","--outname",dest="output",type="str",
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    if not inputfile or not output:
        optparser.print_help()
        sys.exit(1)
    
    fraglen(inputfile,output)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


