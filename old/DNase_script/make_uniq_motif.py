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
def uniq_motif(inputmotif,outputmotif,column):
    inf = open(inputmotif)
    outf = open(outputmotif,'w')
    last = inf.readline().strip().split("\t")
    for line in inf:
        ll = line.strip().split("\t")
        if ll[0] == last[0] and int(ll[1])<int(last[2]):
            if float(ll[column-1]) < float(last[column-1]):
                pass
            else:
                last = ll    
        else:
            outf.write("\t".join(last)+"\n")
            last = ll
    outf.write("\t".join(last)+"\n")
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
    optparser.add_option("-i","--inputmotif",dest="inputmotif",type="str",
                         help="input motif bed file , motif file should be sorted by genomic location")
    optparser.add_option("-o","--outputmotif",dest="outputmotif",type="str",
                         help="output motif bed file after uniq")
    optparser.add_option("-c","--column",dest="column",type="int",default = 5,
                         help="column for rank input motif, when 2 region overap , will remove the one with smaller score")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputmotif = options.inputmotif
    outputmotif = options.outputmotif 
    column = options.column
    
    if not inputmotif:
        optparser.print_help()
        sys.exit(1)
        
    uniq_motif(inputmotif,outputmotif,column)
    

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


