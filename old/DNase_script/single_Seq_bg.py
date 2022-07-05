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
def read_bg_matrix(matrixfile):
    inf = open(matrixfile)
    dall_p = {}
    dall_n = {}
    for line in inf:
        ll = line.split()
        dall_p[ll[0]]=ll[1]
        dall_n[ll[0]]=ll[2]
    inf.close()
    return dall_p,dall_n
    
def single_scan(inputfile,outputfile,bgfile):
    dall_p,dall_n = read_bg_matrix(bgfile)
    inf = open(inputfile)
    outf = open(outputfile,'w')
    for line in inf:
        seq = line.strip()
        pv = []
        nv = []
        for i in range(len(seq)-6-1):
            nbg = dall_n[seq[i:(i+6)]]
            pbg = dall_p[seq[(i+1):(i+1+6)]]
            pv.append(pbg)
            nv.append(nbg)
            
        newll = [seq] + pv + nv
        outf.write("\t".join(newll)+"\n")
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
    optparser.add_option("-b","--bg",dest="bg",type="str",
                         help="")
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outputfile",dest="out",type="str",
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    bg = options.bg
    out = options.out
    inputfile = options.inputfile
    if not bg:
        optparser.print_help()
        sys.exit(1)
    
    single_scan(inputfile,out,bg)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

