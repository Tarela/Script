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

def generate_reads_fa(genomefa,output,seqlen):
    chromSeq = ""
    genome = {}
    inf  = open(genomefa)
    outf = open(output,'w')
    for line in inf:
        if line.startswith(">"):
            if not chromSeq == "":
                genome[chrm] = chromSeq
            chrm = line.split()[0].strip(">")
            chromSeq = ""
        else:
            chromSeq += line.strip()
    genome[chrm] = chromSeq
    inf.close()
    count = 0
    for chrom in genome:
        sequence = genome[chrom]
        for i in range(len(sequence)-seqlen):
            
            read = sequence[i:(i+seqlen)].upper()
            if "N" in read:
                continue
            count += 1
            outf.write(">r"+str(count)+"\n"+read+"\n")
    outf.close()
        #print chrom,len(genome[chrom])


            
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
    optparser.add_option("-l","--readlen",dest="readlen",type="int",
                         help="")
 #   optparser.add_option("-n","--outn",dest="outn",type="str",
 #                        help="")
                       
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    out = options.out
  #  outn = options.outn

    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    generate_reads_fa(inputfile,out,options.readlen)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



