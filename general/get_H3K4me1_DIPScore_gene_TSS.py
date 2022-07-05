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
import numpy
#from CistromeAP.taolib.CoreLib.BasicStat.Func import *
#from CistromeAP.jianlib.BwReader import BwIO



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
def bwsig_Ave(bwfile,chrm,start,end,points):
    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,1)
#    print sp(cmd)
    bwS = sp(cmd)[0].strip()
    if bwS == "":
        return 'na'
    else:
        sigAve = float( bwS)#numpy.array(map(float,sp(cmd)[0].strip().split("\t")))
        return sigAve#[CpGcount,aveME]

def DIPscore(BWH,chrm,TSS,strand):
    ## by default: left = -5000 ~ -500
    ##             center = -500 ~ 1500
    ##             right = 1500 ~ 5000
    if strand == "+":
        upstream = bwsig_Ave(BWH,chrm,max(0,TSS-5000),max(0,TSS-500),1)
        Center = bwsig_Ave(BWH,chrm,max(0,TSS-500),TSS+1500,1)
        downstream = bwsig_Ave(BWH,chrm,max(0,TSS+1500),TSS+5000,1)
    else:
        upstream = bwsig_Ave(BWH,chrm,max(0,TSS+500),TSS+5000,1)
        Center = bwsig_Ave(BWH,chrm,max(0,TSS-1500),TSS+500,1)
        downstream = bwsig_Ave(BWH,chrm,max(0,TSS-5000),max(0,TSS-1500),1)

    if upstream == "na" or Center == "na" or downstream == "na":
        DIP = "na"
    else:
        DIP = -1*numpy.log10((Center+0.001) / ((upstream + downstream+0.001)/4))

    return DIP


def get_signal(inputfile,output,signalbw,bwfolder):
    signalbw = signalbw.strip().strip(',').split(',')
    if not bwfolder:
        bwfolder = ""
    if not bwfolder.endswith('/'):
        bwfolder += '/'
    bwHs = []
    for sb in signalbw:
        if "/" in sb:#sb.startswith('/mnt/Storage'):
            bwHs.append(sb)
        else:
            bwHs.append(bwfolder + sb)
    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        if "_" in line.split()[0]:
            continue
        ll = line.split()
        if ll[5] == "+":
            TSS = int(ll[1])
        else:
            TSS = int(ll[2])
      #  center = ( int(ll[1])+int(ll[2]) ) /2
        
        for bwHandle in bwHs:
            addsig = DIPscore(bwHandle,ll[0],TSS,ll[5])
            #addsig  = bwsig_Ave(bwHandle,ll[0],max(0,center-extend),center+extend,1)
            ll.append(addsig)
        outf.write("\t".join(map(str,ll))+"\n")
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
    optparser.add_option("-i","--input",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
#    optparser.add_option("-w","--signalbw",dest="bw",type="str",default = "ATAC_64.bw,ATAC_256.bw,ATAC_1k.bw,ATAC_oblong.bw,ATAC_dome.bw,ATAC_24hpf.bw,DNase_1cell.bw,DNase_256.bw,DNase_1k.bw,DNase_Dome.bw",
#                         help="")
    optparser.add_option("-w","--signalbw",dest="bw",type="str",
                         help="")

    optparser.add_option("--bwfolder",dest="bwfolder",type="str",default=".",
                         help="")
#    optparser.add_option("--strand",dest="bwfolder",type="str",default=0,
#                         help="consider strand info(6th col) or not, 1 : consider strand, 0 : not consider, default=0")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    signalbw = options.bw
    bwfolder = options.bwfolder
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,signalbw,bwfolder)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


