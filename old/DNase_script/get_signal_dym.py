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
from CistromeAP.taolib.CoreLib.BasicStat.Func import *
from CistromeAP.jianlib.BwReader import BwIO
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

def get_signal(inputfile,output,vp,vm,dp,dm):
    p=BwIO(vp)
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    vpBw=BigWigFile(open(vp, 'rb'))
    vmBw=BigWigFile(open(vm, 'rb'))
    dpBw=BigWigFile(open(dp, 'rb'))
    dmBw=BigWigFile(open(dm, 'rb'))
    inf = open(inputfile)
    outf = open(output,'w')
    colnames = ["chrom","start","end","seq","motifscore","strand","LncapARsignal","LncapDNaseCutsite","LncapDNaseFrag","K562DNaseFrag","LncapFP","K562FP","overARpeak","VehPlus","VehMinus","DHTPlus","DHTMinus"]
    outf.write("\t".join(colnames)+"\n")
    for line in inf:
        if line.startswith("chrom"):
            continue
        ll = line.split()
        if not chrom_len.has_key(ll[0]):
            continue
        signal=vpBw.summarize(ll[0],int(ll[1])-50,int(ll[2])+50,1)
        ll.append(str(float(signal.sum_data)))
        signal=vmBw.summarize(ll[0],int(ll[1])-50,int(ll[2])+50,1)
        ll.append(str(float(signal.sum_data)))
        signal=dpBw.summarize(ll[0],int(ll[1])-50,int(ll[2])+50,1)
        ll.append(str(float(signal.sum_data)))
        signal=dmBw.summarize(ll[0],int(ll[1])-50,int(ll[2])+50,1)
        ll.append(str(float(signal.sum_data)))
        outf.write("\t".join(ll)+"\n")
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
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")
    optparser.add_option("--vehplus",dest="vp",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/Lncap_DNase/Veh_DHT/Raw/Lncap_Veh_cut_plus_normalize.bw",
                         help="")
    optparser.add_option("--vehminus",dest="vm",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/Lncap_DNase/Veh_DHT/Raw/Lncap_Veh_cut_minus_normalize.bw",
                         help="")
    optparser.add_option("--dhtplus",dest="dp",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/Lncap_DNase/Veh_DHT/Raw/Lncap_DHT_cut_plus_normalize.bw",
                         help="")
    optparser.add_option("--dhtminus",dest="dm",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/Lncap_DNase/Veh_DHT/Raw/Lncap_DHT_cut_minus_normalize.bw",
                         help="")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    vp = options.vp
    vm = options.vm
    dp = options.dp
    dm = options.dm
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,vp,vm,dp,dm)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

