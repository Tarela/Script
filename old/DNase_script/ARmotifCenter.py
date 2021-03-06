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
def reverse(seq):
    outseq = []
    for bp in seq:
        if bp.upper() == 'A':
            outseq.append('T')
        elif bp.upper() == 'T':
            outseq.append('A')
        elif bp.upper() == 'C':
            outseq.append('G')
        elif bp.upper() == 'G':
            outseq.append('C')
        elif bp.upper() == 'N': 
            outseq.append('N')
        else :
            print seq,bp
            sys.exit(1)
    return "".join(outseq[::-1])
def ARcenter(inputfile,outfile):
    typedict = {}
    inf = open(inputfile)
    for line in inf:
        ll = line.split()
        if ll[5] == '+':
            seq = ll[3].upper()
        else:
            seq = reverse(ll[3])
        centerseq = seq[7:10]
        if not typedict.has_key(centerseq):
            typedict[centerseq] = 0
        typedict[centerseq] += 1
    inf.close()
    #tt = iteritem(typedict)
    outf = open(outfile,'w')
    for i in sorted(typedict.items(),key=lambda x:x[1],reverse=True):
        outf.write("\t".join(map(str,i))+"\n")
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
                         help="ARmotif file,6+  column")
    optparser.add_option("-o","--outfile",dest="outfile",type="str",
                         help="outfile , 3bp type + count")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    ARcenter(inputfile,options.outfile)
    

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


