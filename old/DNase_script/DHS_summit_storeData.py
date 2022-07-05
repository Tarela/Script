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
import math,time
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

def sitepro_scan(peak,out,w_plus,w_minus):

    inf = open(peak)
    p=BwIO(w_plus)
    q=BwIO(w_minus)
    chrom_len1 = {}
    chrom_len2 = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len1[i['key']] = i['chromSize']
    for i in q.chromosomeTree['nodes']:
        chrom_len2[i['key']] = i['chromSize']
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    footprint = []
    outf = open(out,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        try:
            START = int(ll[1])
        except:
            print 'start:',ll[1],line
            continue
        try:
            END = int(ll[2])
        except:
            print 'end:',ll[2],line
            continue
        if chrom_len1.has_key(ll[0])  and chrom_len2.has_key(ll[0]):
            try:
                DNase100p = sum(list(w_plus_H.summarize(ll[0],(START+END)/2-50,(START+END)/2+50,2).sum_data))
            except:
                #print ll[0],START,END,(START+END)/2-50,(START+END)/2+50,type(START),type(END),line
                DNase100p = sum(list(w_plus_H.summarize(ll[0],(START+END)/2-50,(START+END)/2+50,1).sum_data))
            try:
                DNase100m = sum(list(w_minus_H.summarize(ll[0],(START+END)/2-50,(START+END)/2+50,2).sum_data))
            except:
                DNase100m = sum(list(w_minus_H.summarize(ll[0],(START+END)/2-50,(START+END)/2+50,1).sum_data))
            DNase200p = sum(list(w_plus_H.summarize(ll[0],(START+END)/2-100,(START+END)/2+100,2).sum_data))
            DNase200m = sum(list(w_minus_H.summarize(ll[0],(START+END)/2-100,(START+END)/2+100,2).sum_data))
            DNase300p = sum(list(w_plus_H.summarize(ll[0],(START+END)/2-150,(START+END)/2+150,2).sum_data))
            DNase300m = sum(list(w_minus_H.summarize(ll[0],(START+END)/2-150,(START+END)/2+150,2).sum_data))
            DNase400p = sum(list(w_plus_H.summarize(ll[0],(START+END)/2-200,(START+END)/2+200,2).sum_data))
            DNase400m = sum(list(w_minus_H.summarize(ll[0],(START+END)/2-200,(START+END)/2+200,2).sum_data))
          #  Chip = float(w_chip_H.summarize(ll[0],START,END,1).sum_data)
            p_sum = list(w_plus_H.summarize(ll[0],START-200,START+200,400).sum_data)
            m_sum = list(w_minus_H.summarize(ll[0],START-200,START+200,400).sum_data)
            fpline = (ll+[DNase100p,DNase100m,DNase200p,DNase200m,DNase300p,DNase300m,DNase400p,DNase400m]+p_sum+m_sum) 
            newline = "\t".join(map(str,fpline))+"\n"
            outf.write(newline)
    outf.close()
   # outf = open(out,'w')
   # for fp in footprint:
   #     newline = "\t".join(map(str,fp))+"\n"
   #     outf.write(newline)
   # outf.close()
              #  print score

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",
                         help="")
#========minor options=============


    (options,args) = optparser.parse_args()

    interval = options.interval
    output = options.output
    w_plus = options.w_plus
    w_minus = options.w_minus
    
    sitepro_scan(interval,output,w_plus,w_minus)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

