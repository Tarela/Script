#!/usr/bin/env python2.7
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
import time
#import numpy
#import twobitreader
import gzip
from copy import deepcopy
import Levenshtein 

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------


def dupUMI(umi_list):
    outumi = []
    while 1:
        umi = umi_list[0]
        outumi.append(umi)
        umi_list.remove(umi)
        copy_umi_list  = deepcopy(umi_list)
        for cmp_umi in copy_umi_list:
            if Levenshtein.distance(umi, cmp_umi) <= 1:
                umi_list.remove(cmp_umi)
        if len(umi_list)==0:
            break
    return outumi

def exonOV_rpkm(infile,outname,exonfile,exonlen):

    os.system("bamToBed -i %s -split > %s_split.bed"%(infile,outname))
    os.system("intersectBed -a %s_split.bed -b %s -wao -s > %s_ovexon.bed"%(outname,exonfile,outname))

    read2gene = {}
    inf = open("%s_ovexon.bed"%(outname))
    for line in inf:
        ll = line.split()
        Rname = ll[3]
        Gnames = ll[10].split(",")
        if Gnames == ["-1"]  :
            continue

        if not read2gene.has_key(Rname):
            read2gene[Rname] = []
        for Gname in Gnames:
            if not Gname in read2gene[Rname]:
                read2gene[Rname].append(Gname)

    inf.close()

    gene2UMI = {}
    for read in read2gene.keys():
        if len(read2gene[read]) == 1:
#            print read, read2gene[read]#continue
            gene = read2gene[read][0]
            if not gene2UMI.has_key(gene):
                gene2UMI[gene] = []
            gene2UMI[gene].append(read.split(":")[-1])

    inf = open(exonlen)
    outf = open("%s_counts.txt"%(outname),'w')
    for line in inf:
        ll = line.split()
        if gene2UMI.has_key(ll[0]):
            allUMI = sorted(gene2UMI[ll[0]])
            allUMIcount = len(allUMI)
            uniqUMIcount = len(dupUMI(allUMI))
        else:
            allUMIcount = 0
            uniqUMIcount = 0
        newll = ll+ [ allUMIcount, uniqUMIcount]
        outf.write("\t".join(map(str,newll))+"\n")
         #   print read,read2gene[read]
          #  Nread_multigene += 1
    inf.close()
    outf.close()
    #print Nread,Nread_multigene


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============


    optparser.add_option("-i",dest="inputfile",type="str",
                         help="inputfile, readsXexon ov file")
    optparser.add_option("-o",dest="outname",type="str",
                         help="name of outputfile. the final count table: 4column, geneID;exonlen;totalUMI;uniqUMI")

#    optparser.add_option("--p5",dest="p5",type="str",
#                         help="list of P5barcodes")
#    optparser.add_option("--p7",dest="p7",type="str",
#                         help="list of P7barcodes")
    optparser.add_option("-e",dest="exonfile",type="str",default="/scratch/sh8tv/Data/refgenes/exon_annotation/mm10/mm10_exon_geneID.bed",
                         help="exon annotation file, 6col: chrm;start;end;exonID;relatedgene;strand")
    optparser.add_option("-l",dest="exonlen",type="str",default="/scratch/sh8tv/Data/refgenes/exon_annotation/mm10/mm10_geneID_exonlength.txt",
                         help="exonlength file, 2column: geneID;exonlen")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    exonOV_rpkm(inputfile,options.outname,options.exonfile,options.exonlen)#,options.p5,options.p7,options.startpos)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)




#exonOV_rpkm("mESC_b9sc14_newovexp_raw.bed","testout.bed","/scratch/sh8tv/Data/refgenes/exon_annotation/mm10/mm10_geneID_exonlength.txt")

























