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
import time,scipy
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def caculate_footprint(digital_list,central_max,central_min,flanking_max,flanking_min,cutoff):
    footprint = []
    last_start = "NA"
    last_end = "NA"
    last_value = "NA"
    for i in range(len(digital_list)): ### i is current position -> right bound of central region
        if i < flanking_min or len(digital_list) - flanking_min - central_min < i :
            continue
        for central_length in range(central_min,central_max):
            for flanking_length in range(flanking_min,flanking_max):
                l_end = i
                c_start = i
                if len(digital_list) - flanking_min > i + central_length:
                    c_end = i + central_length
                else:
                    c_end = len(digital_list) - flanking_min
                r_start = c_end
                r_end = r_start + flanking_length
                if r_end > len(digital_list):
                    r_end = len(digital_list)
                    l_start = i - ( r_end - r_start )
                elif i <= flanking_length:
                    l_start = 0
                    r_end = r_start + i
                else:
                    l_start = i - flanking_length

                L = scipy.mean(digital_list[l_start:l_end])
                C = scipy.mean(digital_list[c_start:c_end])
                R = scipy.mean(digital_list[r_start:r_end])
                if L <= C or R <= C:
                    continue
                else:
                    FOS = (C+1)/L + (C+1)/R
                if last_start == "NA" :
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
                elif c_start >= last_end:
                    if last_value < cutoff :
                        footprint.append([last_start,last_end,last_value])
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
                elif FOS < last_value:
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
    if last_value < cutoff:
        footprint.append([last_start,last_end,last_value])
    return footprint

def summary(bwfile,bedfile,out,central_max,central_min,flanking_max,flanking_min,cutoff):
    total_result = []
    p=BwIO(bwfile)
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
    bwHandle=BigWigFile(open(bwfile, 'rb'))
    inf = open(bedfile)
    outf = open(out,'w')
    t = time.time()
    for line in inf:
        ll = line.split()
        if chrom_len.has_key(ll[0]):
            #t = time.time()
            summary = bwHandle.summarize(ll[0],int(ll[1]),int(ll[2]),(int(ll[2])-int(ll[1])))
    #        print "bw sum time",time.time()-t
     #       t=time.time()
            digital = list(summary.sum_data)
      #      print "trans to list time",time.time()-t
      #      t=time.time()
            FT=(caculate_footprint(digital,central_max,central_min,flanking_max,flanking_min,cutoff))
       #     print "scan footprint time",time.time()-t
       #     time.time()
            for ft in FT:
                bed ="\t".join(map(str,[ll[0],int(ll[1])+ft[0],int(ll[1])+ft[1],ll[3],ft[2]]))+"\n"
                outf.write(bed)
            #print "single time",time.time()-t
            #print (int(ll[2])-int(ll[1]))#*1.0/(time.time()-t)
    inf.close()   
    outf.close()
    print "scaning 1st ",time.time()-t

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
                         help="")
#    optparser.add_option("-x","--bwfile_add",dest="bwfile_add",type="str",
#                         help="for multi input ,sep by comma")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")
#    optparser.add_option("--topnumber",dest="topnumber",type="int",default=10000,
#                         help="")
#========minor options============
    optparser.add_option("--c_max",dest="c_max",type="int",default = 40,
                         help="max length of central region,defult 40")
    optparser.add_option("--c_min",dest="c_min",type="int",default = 8,
                         help="min length of central region,default 8")
    optparser.add_option("--f_max",dest="f_max",type="int",default = 6,
                         help="max length of flanking region,default 6")
    optparser.add_option("--f_min",dest="f_min",type="int",default = 5,
                         help="min length of flanking region,default 5")
    optparser.add_option("-c","--cutoff",dest="cutoff",type="float",default = 1,
                         help="cutoff of the footprint score, default 1")

    (options,args) = optparser.parse_args()

    cutoff = options.cutoff
    c_max = options.c_max
    c_min = options.c_min
    f_max = options.f_max
    f_min = options.f_min

    if not options.bedfile or not options.bwfile:# or not bwfile2:
        optparser.print_help()
        sys.exit(1)

    bedfile = options.bedfile
    bwfile = options.bwfile
    out = options.out
    if not bedfile or not bwfile:
        optparser.print_help()
        sys.exit(1)
    summary(bwfile,bedfile,out,c_max,c_min,f_max,f_min,cutoff)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

