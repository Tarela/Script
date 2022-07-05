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

def fraglen(inputfile,output,uplim):
    inf = open(inputfile)
    D = {}
    for i in range(int(float(uplim))):
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
                         help="input file in more than 4 column bed file , format should be [chrm,fragment_start,fragment_end,fragment_length], fragment refer to the combination of 5'end tag and 3'end tag, which is pairedend  fragment")
    optparser.add_option("-o","--outname",dest="output",type="str",default='NA',
                         help="name of your output file, you will have 2 output file named as youname.r and yourname.pdf")

#========minor options=============

    optparser.add_option("--upperlimit",dest="uplim",type="int",default =600,
                         help="upper limit of your paired end size range, default is 600")

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    
    if not inputfile or not output:
        optparser.print_help()
        sys.exit(1)
    
    fraglen(inputfile,output,options.uplim)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


