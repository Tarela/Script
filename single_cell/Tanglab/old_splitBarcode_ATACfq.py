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
import gzip
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def splitATACfq(inputname,outname,folder,p5,p7,startpos,barcodeLen):

    fq1 = "%s/%s/%s_R1.fastq.gz"%(folder,inputname,inputname)
    fq2 = "%s/%s/%s_R2.fastq.gz"%(folder,inputname,inputname)
    inf1 = gzip.open(fq1,"rb")
    inf2 = gzip.open(fq2,"rb")

    line_count = 0
    reads_count = {}

    outf1 = open("%s_R1.fastq"%(outname),"w")
    outf2 = open("%s_R2.fastq"%(outname),"w")

    while 1:
        line_count += 1
        p1 = inf1.readline()
        p2 = inf2.readline()

        if line_count%4 == 1:
            p1l1 = p1
            p2l1 = p2
        elif line_count%4 == 2:
            p1l2 = p1[startpos:]
            p2l2 = p2[startpos:]
            line_p5 = p1[:barcodeLen]
            line_p7 = p2[:barcodeLen] 
            #print "p1",p1,line_p5
            #print "p2",p2,line_p7               
        elif line_count%4 == 3:
            p1l3 = p1
            p2l3 = p2
        else:
            p1l4 = p1[startpos:]
            p2l4 = p2[startpos:]

            barcode_key = line_p5 + "_" + line_p7
            if not reads_count.has_key(barcode_key):
                reads_count[barcode_key] = 0
            reads_count[barcode_key] += 1

            if line_p5 == p5 and line_p7 == p7:
                outf1.write(p1l1)
                outf1.write(p1l2)
                outf1.write(p1l3)
                outf1.write(p1l4)
                outf2.write(p2l1)
                outf2.write(p2l2)
                outf2.write(p2l3)
                outf2.write(p2l4)

        if p1.strip() == "" and p2.strip() == "":
            break
    inf1.close()
    inf2.close()
    outf1.close()
    outf2.close()

    outf = open("%s_BCcount.txt"%outname,'w')
    for BC,count in sorted(reads_count.iteritems(),key=lambda x:x[1],reverse=True):
        outf.write("%s\t%s\n"%(BC,count))
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


    optparser.add_option("-i",dest="inputname",type="str",
                         help="")
    optparser.add_option("-o",dest="outname",type="str",
                         help="")
    optparser.add_option("--p5",dest="p5",type="str",
                         help="")
    optparser.add_option("--p7",dest="p7",type="str",
                         help="")
    optparser.add_option("--startpos",dest="startpos",type="int",default=27,
                         help="")
    optparser.add_option("--barcodeLen",dest="barcodeLen",type="int",default=8,
                         help="")
    optparser.add_option("-f","--folder",dest="folder",type="str",default="/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/05_20181030_mESC_tagmentation_CAbeads/combined_ATAC/",
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputname = options.inputname
    if not inputname:
        optparser.print_help()
        sys.exit(1)
    splitATACfq(inputname,options.outname,options.folder,options.p5,options.p7,options.startpos,options.barcodeLen)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


