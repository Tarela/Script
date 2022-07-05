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

def scan_fp(plusdnase,minusdnase,bed,out,upstream,downstream):
    p=BwIO(plusdnase)
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    bwHandle1=BigWigFile(open(plusdnase, 'rb'))
    bwHandle2=BigWigFile(open(minusdnase, 'rb'))
    inf = open(bed)
    outf = open(out,'w')
    for line in inf:
        ll = line.split()
        if not chrom_len.has_key(ll[0]):
            continue
        if int(ll[1])<upstream : 
            continue
        signal1=bwHandle1.summarize(ll[0],int(ll[1])-upstream,int(ll[2])+downstream,(int(ll[2])+downstream-int(ll[1])+upstream))
        signal2=bwHandle2.summarize(ll[0],int(ll[1])-upstream,int(ll[2])+downstream,(int(ll[2])+downstream-int(ll[1])+upstream))
        #ll.append(str(float(signal.sum_data)))
        newll = ll[:6]+map(str,list(signal1.sum_data))+map(str,list(signal2.sum_data))
        outf.write("\t".join(newll)+"\n")
    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """for scan footprint on given motif region , current only support plus strand motif"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-p","--plusdnase",dest="plusdnase",type="str",
                         help="")
    optparser.add_option("-m","--minusdnase",dest="minusdnase",type="str",
                         help="")
    optparser.add_option("-b","--bed",dest="bed",type="str",
                         help="current only plus strand motif")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")
#========minor options=============
    optparser.add_option("--upstream",dest="upstream",type="int",default=50,
                         help="up means up from motif start")
    optparser.add_option("--downstream",dest="downstream",type="int",default=50,
                         help="down from motif end")


    (options,args) = optparser.parse_args()

    plusdnase = options.plusdnase
    minusdnase = options.minusdnase
    bed = options.bed
    out = options.out
    upstream = options.upstream
    downstream = options.downstream
   
    if not plusdnase or not minusdnase or not bed or not out:
        optparser.print_help()
        sys.exit(1)
    
    scan_fp(plusdnase,minusdnase,bed,out,upstream,downstream)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

