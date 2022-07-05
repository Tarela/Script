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

    
def getsignal(inputfile,outputfile,pcut,DHT,Veh,pspan):

    
#    p=BwIO(pcut)
#    chrom_len = {}
#    for i in p.chromosomeTree['nodes']:
#        chrom_len[i['key']] = i['chromSize']
    pcutbw = BigWigFile(open(pcut, 'rb'))
    dht = BigWigFile(open(DHT, 'rb'))
    veh = BigWigFile(open(Veh, 'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    inf.seek(0)
    outf = open(outputfile,'w')

    for line in inf:
        ll = line.split()
#        if not chrom_len.has_key(ll[0]):
#            continue
        cut = list(pcutbw.summarize(ll[0],int(ll[1]) + ml/2 -pspan ,int(ll[1]) + ml/2 +pspan ,2*pspan).sum_data)
        TC = sum(cut)
        C = sum(cut[(pspan-ml/2) : (pspan-ml/2+ml)])
        L = sum(cut[(pspan-ml/2-ml):(pspan-ml/2)])
        R = sum(cut[(pspan-ml/2+ml):(pspan-ml/2+2*ml)])
        FOS = -1*( (C+1)/(R+1) + (C+1)/(L+1) )
        dhtnum = sum(list(dht.summarize(ll[0],int(ll[1]) + ml/2 -pspan ,int(ll[1]) + ml/2 +pspan ,2).sum_data)) + 1
        vehnum = sum(list(veh.summarize(ll[0],int(ll[1]) + ml/2 -pspan ,int(ll[1]) + ml/2 +pspan ,2).sum_data)) + 1
        newll = ll + [TC,FOS,dhtnum,vehnum]
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

    optparser.add_option("-w","--cutbw",dest="cut",type="str",
                         help="")
    optparser.add_option("-d","--DHTbw",dest="dht",type="str",
                         help="")
    optparser.add_option("-v","--Vehbw",dest="veh",type="str",
                         help="")


    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 100,
                         help="window size for tagcount and dymDHS, default = 100 , means total length = 200")

 
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    cutbw = options.cut
    dht = options.dht
    veh = options.veh
    pspan = options.pspan
    
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,cutbw,dht,veh,pspan)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
