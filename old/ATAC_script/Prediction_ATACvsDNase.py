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

    
def getsignal(inputfile,outputfile,ATAC100,ATAC247,ATACall,DNase,pspan):

    
#    p=BwIO(pcut)
#    chrom_len = {}
#    for i in p.chromosomeTree['nodes']:
#        chrom_len[i['key']] = i['chromSize']
    ATAC100bw = BigWigFile(open(ATAC100, 'rb'))
    ATAC247bw = BigWigFile(open(ATAC247, 'rb'))
    ATACallbw = BigWigFile(open(ATACall, 'rb'))
    DNasebw = BigWigFile(open(DNase, 'rb'))

    inf = open(inputfile)    
    outf = open(outputfile,'w')

    for line in inf:
        ll = line.split()
        if ll[0] == 'chrY':
            continue
#        print [ll[0],(int(ll[1])+int(ll[2]))/2 -pspan ,(int(ll[1])+int(ll[2]))/2 -pspan]
        ATAC100_signal = float(ATAC100bw.summarize(ll[0],(int(ll[1])+int(ll[2]))/2 -pspan ,(int(ll[1])+int(ll[2]))/2 +pspan,1).sum_data)/(2*pspan)
        ATAC247_signal = float(ATAC247bw.summarize(ll[0],(int(ll[1])+int(ll[2]))/2 -pspan ,(int(ll[1])+int(ll[2]))/2 +pspan,1).sum_data)/(2*pspan)
        ATACall_signal = float(ATACallbw.summarize(ll[0],(int(ll[1])+int(ll[2]))/2 -pspan ,(int(ll[1])+int(ll[2]))/2 +pspan,1).sum_data)/(2*pspan)
        DNase_signal = float(DNasebw.summarize(ll[0],(int(ll[1])+int(ll[2]))/2 -pspan ,(int(ll[1])+int(ll[2]))/2 +pspan,1).sum_data)/(2*pspan)
        

        newll = ll + [ATAC100_signal,ATAC247_signal,ATACall_signal,DNase_signal]
        outf.write("\t".join(map(str,newll))+"\n")

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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")                         

    optparser.add_option("-a","--atac100bw",dest="atac100",type="str",
                         help="")
    optparser.add_option("-b","--atac247bw",dest="atac247",type="str",
                         help="")
    optparser.add_option("-c","--atacallbw",dest="atacall",type="str",
                         help="")
    optparser.add_option("-d","--dnasebw",dest="dnase",type="str",
                         help="")

    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 100,
                         help="window size for tagcount and dymDHS, default = 100 , means total length = 200")

 
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    atac100 = options.atac100
    atac247 = options.atac247
    atacall = options.atacall
    dnase = options.dnase
    pspan = options.pspan
    
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,atac100,atac247,atacall,dnase,pspan)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
