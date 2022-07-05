'''
Created on 2011-5-24

@author: undergraduate
'''
#!/usr/bin/env python
#Time-stamp:<2011-05-24 Sheng'en>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import logging
from optparse import OptionParser
import subprocess
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

def sudocount(x):
    if x==0:
        return 0.001
    else:
        return x

def Readbw(bwfile,chrm,start,end,n):
    bwHandle=BigWigFile(open(bwfile, 'rb'))
    summary = bwHandle.summarize(chrm,int(start),int(end),(int(end)-int(start))/n)
    count = map(sudocount,summary.valid_count)
    sum = summary.sum_data
    scores = list(sum/count)
    return scores

def summary(bwfile,bedfile,topnumber,out):
    total_result = []
    p=BwIO(bwfile)
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    bwHandle=BigWigFile(open(bwfile, 'rb'))
    inf = open(bedfile)
    t = time.time()
    for line in inf:
        ll = line.split()
        ll[3]="-"
        if chrom_len.has_key(ll[0]):
            summary = bwHandle.summarize(ll[0],int(ll[1]),int(ll[2]),1)
            if summary.valid_count == 0:
                mean_value = 0
            else:
                mean_value = (summary.sum_data/summary.valid_count)[0]
            total_result.append(ll+[mean_value])
    inf.close()   
    total_result.sort(reverse=True,key=lambda x:x[-1])
    outf = open(out,'w')
    print "scaning 1st ",time.time()-t
    t=time.time()
    for i in range(topnumber):
        ll = total_result[i]
        summary = bwHandle.summarize(ll[0],int(ll[1]),int(ll[2]),(int(ll[2])-int(ll[1])))
        additional_value = ",".join(map(str,list(summary.sum_data)))
        result = map(str,(ll+[additional_value]))
        outf.write("\t".join(result)+"\n")
    outf.close()
    print "scaning 2nd ",time.time()-t
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-p deletepath> "
    description = """delete all unuse pipeline result in given path"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-b","--bedfile",dest="bedfile",type="str",
                         help="")
    optparser.add_option("-w","--bwfile",dest="bwfile",type="str",
                         help="")
    optparser.add_option("-o","--outfile",dest="outfile",type="str",
                         help="")
    optparser.add_option("--topnumber",dest="topnumber",type="int",default=10000,
                         help="")
    (options,args) = optparser.parse_args()

    bedfile = options.bedfile
    bwfile = options.bwfile
    outfile = options.outfile
    topnumber = options.topnumber
    if not bedfile or not bwfile:
        optparser.print_help()
        sys.exit(1)
    
    summary(bwfile,bedfile,topnumber,outfile)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

