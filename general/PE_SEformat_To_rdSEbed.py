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
import string
import random

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
def rdSE(inputfile,outputfile):
    random.seed(1)
    inf = open(inputfile)
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()
        a = random.randint(0,1)
        if len(ll) >= 6:
            V1 = ll[3]
            V2 = ll[4]
        else:
            V1 = "."
            V2 = "."
        if a == 0:
            newll = [ll[0],int(float(ll[1])),int(float(ll[1]))+1,V1,V2,"+"] 
        else:
            newll = [ll[0],int(float(ll[2]))-1,int(float(ll[2])),V1,V2,"-"] 
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

    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    if not inputfile:
        optparser.print_help()
        sys.exit(1)

    rdSE(inputfile,output)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

