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
import math
from cistrome import regions as c
import twobitreader
#from CistromeAP.taolib.CoreLib.BasicStat.Func import *
#from CistromeAP.jianlib.BwReader import BwIO
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

    
def make_chip(chipbw_handle,ll,span):
    center = (int(ll[1])+int(ll[2]))/2
    total = (chipbw_handle.summarize(ll[0],max(0,center-span),center+span,2*span/10).sum_data)
    vd = (chipbw_handle.summarize(ll[0],max(0,center-span),center+span,2*span/10).valid_count)
    outraw = total/(vd+0.0001)
    out = list(outraw)
    return out
    

    
def getsignal(inputfile,outputfile,dnasebw,atacbw1,atacbw2,atacbw3,conserve,chipbw,dnase_span):

    
 #   p=BwIO(pcut)
 #   chrom_len = {}
 #   for i in p.chromosomeTree['nodes']:
 #       chrom_len[i['key']] = i['chromSize']
    conservebw = BigWigFile(open(conserve, 'rb'))
    chip_bw = BigWigFile(open(chipbw,'rb'))
    dnase_bw = BigWigFile(open(dnasebw,'rb'))
    atac_dhs = BigWigFile(open(atacbw1,'rb'))
    atac_mono_dhs = BigWigFile(open(atacbw2,'rb'))
    atac_all = BigWigFile(open(atacbw3,'rb'))
    
    #atacnuc_bw = BigWigFile(open(atacnucbw,'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    #ml = int(testll[2]) - int(testll[1])
    
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()
        if ll[0] == 'chrY':
            continue
        conserveout = make_chip(conservebw,ll,dnase_span)
        chipout = make_chip(chip_bw,ll,dnase_span)
        dnaseout = make_chip(dnase_bw,ll,dnase_span)
        atac_dhs_out = make_chip(atac_dhs,ll,dnase_span)
        atac_mono_dhs_out = make_chip(atac_mono_dhs,ll,dnase_span)
        atac_all_out = make_chip(atac_all,ll,dnase_span)
        
        newll = ll + atac_dhs_out + atac_mono_dhs_out + atac_all_out + dnaseout + chipout + conserveout
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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")                         
    optparser.add_option("-a","--atacbw1",dest="atac_dhs",type="str",
                         help="")
    optparser.add_option("-b","--atacbw2",dest="atac_mono_dhs",type="str",
                         help="")
    optparser.add_option("-c","--atacbw3",dest="atac_all",type="str",
                         help="")
    optparser.add_option("-d","--dnasebw",dest="dnase",type="str",
                         help="")

    optparser.add_option("--conserve",dest="conserve",type="str",default = "/mnt/Storage/home/huse/Data/conserv/hg19_v_conserv.bw",
                         help="default is /mnt/Storage/home/huse/Data/conserv/hg19_v_conserv.bw")
    optparser.add_option("--chipbw",dest="chipbw",type="str",
                         help="")
                         
    optparser.add_option("--Dspan",dest="Dspan",type="int",default = 200,
                         help="range of signal,1bp per bin, default = 200 , means total length = 400")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    dnase_span = options.Dspan
    conserve = options.conserve
    chipbw = options.chipbw
    atacbw1 = options.atac_dhs
    atacbw2 = options.atac_mono_dhs
    atacbw3 = options.atac_all
    dnasebw = options.dnase
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    getsignal(inputfile,outputfile,dnasebw,atacbw1,atacbw2,atacbw3,conserve,chipbw,dnase_span)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
