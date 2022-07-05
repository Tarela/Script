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
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf = open(out,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        p_sum = list(w_plus_H.summarize(ll[0],int(ll[1]),int(ll[2]),int(ll[2])-int(ll[1])).sum_data)
        m_sum = list(w_minus_H.summarize(ll[0],int(ll[1]),int(ll[2]),int(ll[2])-int(ll[1])).sum_data)
        fp=(ll+p_sum+m_sum)
        newline = "\t".join(map(str,fp))+"\n"
        outf.write(newline)
    
    outf.close()
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
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_plus_50_100.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_minus_50_100.bw",
                         help="")
#========minor options=============


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)


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

