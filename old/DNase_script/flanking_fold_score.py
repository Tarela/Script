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
import math
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

def get_signal(inputfile,output,Pbw,Nbw,score_range):
    persudo = 0.2
    p=BwIO(Pbw)
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    PH=BigWigFile(open(Pbw, 'rb'))
    NH=BigWigFile(open(Nbw,'rb'))
    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        if not chrom_len.has_key(ll[0]):
            continue
        motif_len = int(ll[2]) - int(ll[1])
        Psignal=list(PH.summarize(ll[0],max(int(ll[1])-100,0),int(ll[1])+100,200).sum_data)
        Nsignal=list(NH.summarize(ll[0],max(int(ll[1])-100,0),int(ll[1])+100,200).sum_data)
        DNase = sum(Psignal) + sum(Nsignal)
        
        if ll[5] == '+':
            S_up_same = sum(Psignal[(100 - score_range):100])
            S_up_diff = sum(Nsignal[(100 - score_range):100])
            S_down_same = sum(Psignal[(100 + motif_len) : 100 + motif_len + score_range])
            S_down_diff = sum(Nsignal[(100 + motif_len) : 100 + motif_len + score_range])
        
        elif ll[5] == '-':
            S_up_same = sum(Nsignal[(100 + motif_len) : 100 + motif_len + score_range])
            S_up_diff = sum(Psignal[(100 + motif_len) : 100 + motif_len + score_range])
            S_down_same = sum(Nsignal[(100 - score_range):100])
            S_down_diff = sum(Psignal[(100 - score_range):100])
        else:
            print line
            sys.exit(1)
        
    #    if S_up_same == 0 or S_up_diff ==0 or S_down_same == 0 or S_down_diff == 0:
    #        continue
        FPscore1 = math.log(( S_up_same+persudo)*(S_down_diff+persudo) / ((S_up_diff+persudo)*(S_down_same+persudo))  , 2)
        FPscore2 = math.sqrt(S_up_same) + math.sqrt(S_down_diff) - math.sqrt(S_up_diff) - math.sqrt(S_down_same)
        
        
        ll.extend([DNase,FPscore1,FPscore2])
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
    optparser.add_option("-p","--pbw",dest="pbw",type="str",
                         help="")
    optparser.add_option("-n","--nbw",dest="nbw",type="str",
                         help="")
    optparser.add_option("--range",dest="score_range",type="int",default = 20,
                         help="")
                         
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    positive_bw = options.pbw
    negative_bw = options.nbw
    score_range = options.score_range
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,positive_bw,negative_bw,score_range)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

