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
import numpy
import twobitreader
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def judge_same_read(readA,readB):
    llA = readA.split()
    llB = readB.split()
    if llA[0] == llB[0] and llA[1] == llB[1] and llA[2] == llB[2] and llA[5] == llB[5]:
        return 1
    else:
        return 0
        
def filterDup(inputbed, outputfilename, dupN):
    inf = open(inputbed)
    outf = open(outputfilename,'w')
    
    last_reads = "NA"
    now_dup_num = 0
    read_name_idx = 0
    for line in inf:        
        if last_reads == "NA" or judge_same_read(last_reads, line):
            last_reads = line
            now_dup_num += 1
        
        else:
            # write lastreads
            read_name_idx += 1
            write_dup_N = max(min(now_dup_num,dupN),1)
            for n in range(write_dup_N):
                read_name = "r"+str(read_name_idx)+"d"+str(n+1)
                ll = last_reads.split()
                newll = ll[:3]+[read_name,".",ll[5]]
                outf.write("\t".join(map(str,newll))+"\n")

            now_dup_num = 1
            last_reads = line
            #now_dup_num += 1
        
    read_name_idx += 1
    write_dup_N = max(min(now_dup_num,dupN),1)
    for n in range(write_dup_N):
        read_name = "r"+str(read_name_idx)+"d"+str(n+1)
        ll = last_reads.split()
        newll = ll[:3]+[read_name,".",ll[5]]
        outf.write("\t".join(map(str,newll))+"\n")

    outf.close()
    inf.close()



# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="")              
    optparser.add_option("-o","--out",dest="outputfilename",type="str",
                         help="")
    optparser.add_option("-d","--dupN",dest="dupN",type="int",default=1,
                         help="number of kept duplicate, default =1 ")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    if not inputbed:
        optparser.print_help()
        sys.exit(1)
    
    filterDup(inputbed, options.outputfilename, options.dupN)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)




