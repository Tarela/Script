'''
Created on XXXX-XX-XX

@author: Tarela
'''
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

def pipeline(motif,dnasepluscut,dnaseminuscut,dnasepeak,tfpeak,tfsignal,name):
    ### 1. make motif + dnase + plus file
    motif_plus_dnase = name+"_plusS1k_DHS.bed"
    cmd1 = """intersectBed -a %s -b %s -u | awk '{if ($6=="+" && $5 > 1000) print $0;}' > %s"""%(motif,dnasepeak,motif_plus_dnase)
    motif_dnase_peak = name+"_plusS1k_DHS_peak_top5k.bed"
    cmd2 = """intersectBed -a %s -b %s -u | sort -k 5,5gr | head -n 5000 > %s"""%(motif_plus_dnase,tfpeak,motif_dnase_peak)
    motif_plus_dnase_signal = name+"_plusS1k_DHS_TFsignal.bed"
    motif_plus_dnase_signal_top5k = name+"_plusS1k_DHS_noTFsignal_top5k.bed"
    cmd3 = """python /home/sh430/Project/ChongzhiDNase/Script/get_signal.py -i %s -o %s -w %s"""%(motif_plus_dnase,motif_plus_dnase_signal,tfsignal)
    cmd4 = """awk '{if ($7==0) print $0;}' %s | sort -k 5,5gr | head -n 5000 > %s"""%(motif_plus_dnase_signal,motif_plus_dnase_signal_top5k) 
    cmd5 = """siteproBW -w %s -w %s --pf-res=1 --span=50 -l +tag_on_+motif -l -tag_on_+motif -b %s --name=%s"""%(dnasepluscut,dnaseminuscut,motif_dnase_peak,name+"_plusS1k_DHS_peak_top5k_sitepro")
    cmd6 = """siteproBW -w %s -w %s --pf-res=1 --span=50 -l +tag_on_+motif -l -tag_on_+motif -b %s --name=%s"""%(dnasepluscut,dnaseminuscut,motif_plus_dnase_signal_top5k,name+"_plusS1k_DHS_noTFsignal_top5k_sitepro")
    newline = "\n".join([cmd1,cmd2,cmd3,cmd4,cmd5,cmd6])+"\n\n"
    print newline


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-m","--motif",dest="motif",type="str",
                         help="motif bed file")
    optparser.add_option("--dnasepluscut",dest="dnasepluscut",type="str",
                         help="dnase plus cut site , bw file")
    optparser.add_option("--dnaseminuscut",dest="dnaseminuscut",type="str",
                         help="dnase minus cut site , bw file")
    optparser.add_option("--dnasepeak",dest="dnasepeak",type="str",
                         help="macs2 peak from dnase data")
    optparser.add_option("--tfpeak",dest="tfpeak",type="str",
                         help="macs2 peak from tf chipseq data")
    optparser.add_option("--tfsignal",dest="tfsignal",type="str",
                         help="tf signal , bw file")
    optparser.add_option("-n","--name",dest="name",type="str",
                         help="output name")
#========minor options=============

    (options,args) = optparser.parse_args()

    motif = options.motif
    dnasepluscut = options.dnasepluscut
    dnaseminuscut = options.dnaseminuscut
    dnasepeak = options.dnasepeak
    tfpeak = options.tfpeak
    tfsignal = options.tfsignal
    name = options.name
    if not motif:
        optparser.print_help()
        sys.exit(1)
    
    pipeline(motif,dnasepluscut,dnaseminuscut,dnasepeak,tfpeak,tfsignal,name)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

