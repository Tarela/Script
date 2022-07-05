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

def trim_fq(inputname,outLen):
    inf = open(inputname + ".fastq")
    outfile = inputname + "_" + str(outLen) + "bp.fastq"
    outf = open(outfile,'w')
    count = 0
    readlen_last = "NA" 
    countNO = 0
    for line in inf:
        count += 1
        if count % 4 == 1:
            readname = line
            #outf.write(line)
        elif count % 4 == 2:
            readlen = len(line.strip())
            if readlen < outLen:
                WRITE = 0
                countNO += 1
                readseq = line.strip()
            else:
                WRITE = 1
                readseq = line.strip()[:outLen]+"\n"
#                print inputname,readname, line.strip(), readlen, outLen
            #outf.write(line.strip()[:outLen]+"\n")
            #print line.strip(),readlen_this
        elif count % 4 == 3:
            tag = line
        else:
            if WRITE == 1:
                outf.write(readname)
                outf.write(readseq)
                outf.write(tag)
                outf.write(line.strip()[:outLen]+"\n")
            else:
                outf.write(readname)
                outf.write(readseq)
                outf.write(tag)
                outf.write(line)
             
#            if readlen_last == "NA":
#                pass
#            else:
#                if readlen_this == readlen_last:
#                    pass
#                else:
#                    print "in-consistant readlen",name,readlen_this, readlen_last
#            readlen_last = readlen_this
                
#        if count == 200:
#            break
#    print name,readlen_this
    inf.close()
    outf.close()     
    print inputname,countNO*1.0/(count+1),count,countNO,

            
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputname",dest="inputname",type="str",
                         help="")
#    optparser.add_option("-o","--out",dest="out",type="str",
#                         help="")    
    optparser.add_option("-l","--readlen",dest="readlen",type="int",
                         help="")
 #   optparser.add_option("-n","--outn",dest="outn",type="str",
 #                        help="")
                       
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputname = options.inputname
  #  out = options.out
  #  outn = options.outn

    if not inputname:
        optparser.print_help()
        sys.exit(1)
    
    trim_fq(inputname,options.readlen)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



