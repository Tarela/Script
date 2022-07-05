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

def summary(bwfile1,bwfile2,bwfile_add,bedfile,topnumber,out):
    total_result = []
    p=BwIO(bwfile1)
    q=BwIO(bwfile2)
    chrom_len1 = {}
    chrom_len2 = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len1[i['key']] = i['chromSize']
    for i in q.chromosomeTree['nodes']:
        chrom_len2[i['key']] = i['chromSize']
    bwHandle1=BigWigFile(open(bwfile1, 'rb'))
    bwHandle2=BigWigFile(open(bwfile2, 'rb'))
    inf = open(bedfile)
    t = time.time()
    for line in inf:
        ll = line.split()
        ll[3]="-"
        if chrom_len1.has_key(ll[0]) and chrom_len2.has_key(ll[0]):
            summary = bwHandle1.summarize(ll[0],int(ll[1]),int(ll[2]),1)
            if summary.valid_count == 0:
                mean_value1 = 0
            else:
                mean_value1 = (summary.sum_data/summary.valid_count)[0]
            summary = bwHandle2.summarize(ll[0],int(ll[1]),int(ll[2]),1)
            if summary.valid_count == 0:
                mean_value2 = 0
            else:
                mean_value2 = (summary.sum_data/summary.valid_count)[0]
            total_result.append(ll+[mean_value1+mean_value2])
    inf.close()   
    total_result.sort(reverse=True,key=lambda x:x[-1])
    bwHs=[]
    for i in bwfile_add:
        bwHs.append(BigWigFile(open(i, 'rb')))
    outf = open(out,'w')
    print "scaning 1st ",time.time()-t
    t=time.time()
    for i in range(min(len(total_result),topnumber)):
        ll = total_result[i]
        summary = bwHandle1.summarize(ll[0],int(ll[1]),int(ll[2]),(int(ll[2])-int(ll[1])))
        additional_value1 = ",".join(map(str,list(summary.sum_data)))
        summary = bwHandle2.summarize(ll[0],int(ll[1]),int(ll[2]),(int(ll[2])-int(ll[1])))
        additional_value2 = ",".join(map(str,list(summary.sum_data)))
        result = map(str,(ll+[additional_value1,additional_value2]))
        for bwH in bwHs: 
            summary = bwH.summarize(ll[0],int(ll[1]),int(ll[2]),(int(ll[2])-int(ll[1])))
            additional_value_add = ",".join(map(str,list(summary.sum_data)))
            result.append(additional_value_add)
        outf.write("\t".join(result)+"\n")
    outf.close()
    print "scaning 2nd ",time.time()-t

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.<"
    description = """>.<"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-b","--bedfile",dest="bedfile",type="str",
                         help="")
    optparser.add_option("-w","--bwfile",dest="bwfile",type="str",
                         help="input should be plus.bw,minus.bw")
    optparser.add_option("-x","--bwfile_add",dest="bwfile_add",type="str",
                         help="for multi input ,sep by comma")
    optparser.add_option("-o","--outfile",dest="outfile",type="str",
                         help="")
    optparser.add_option("--topnumber",dest="topnumber",type="int",default=10000,
                         help="")
    (options,args) = optparser.parse_args()

    if not options.bedfile or not options.bwfile:# or not bwfile2:
        optparser.print_help()
        sys.exit(1)

    bedfile = options.bedfile
    bwfile1 = options.bwfile.split(",")[0]
    bwfile2 = options.bwfile.split(",")[1]
    bwfile_add = options.bwfile_add.split(",")
    outfile = options.outfile
    topnumber = options.topnumber
    #if not bedfile or not bwfile1 or not bwfile2:
    #    optparser.print_help()
    #    sys.exit(1)
    summary(bwfile1,bwfile2,bwfile_add,bedfile,topnumber,outfile)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

