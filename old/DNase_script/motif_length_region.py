#!/usr/bin/env python
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
import math
import subprocess

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
    
def calculate(motif_dir,outputfile):
    allfiles = os.listdir(motif_dir)
    outf = open(outputfile,'w')
    for f in allfiles:
        if f.endswith("_uniq.bed"):
            name = f.rstrip("_uniq.bed")
            inf = open(motif_dir + "/" + f)
            line = inf.readline()
            motiflen = int(line.split()[2]) - int(line.split()[1])
            motifnum = sp('wc -l '+motif_dir + "/" + f)[0].split()[0]
            newll = [name,motiflen,motifnum]
            outf.write("\t".join(map(str,newll))+"\n")
            inf.close()
    outf.close()
    
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog  "
    description = """"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-d","--inputdir",dest="inputdir",type="str",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    (options,args) = optparser.parse_args()

    inputdir = options.inputdir
    output = options.output
    if not inputdir or not output:
        optparser.print_help()

        sys.exit(1)
    calculate(inputdir,output)



if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

