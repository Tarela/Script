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
import string,numpy
from bx.bbi.bigwig_file import BigWigFile

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

#import subprocess
#def sp(cmd):
#    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
#    ac = a.communicate()
#    return ac
#def bwsig_pattern(bwfile,chrm,start,end,points):
#    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
#    sigPat = sp(cmd)[0].strip().split("\t")
#    return sigPat

def get_signal(inputfile,output,plusBW,minusBW,bwfolder,extend):

    if not bwfolder:
        bwfolder = "./"
    if not bwfolder.endswith('/') :
        bwfolder += '/'

    plus = BigWigFile(open(bwfolder + plusBW, 'rb'))
    minus = BigWigFile(open(bwfolder + minusBW, 'rb'))

    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        if "_" in ll[0]:
            continue
        if len(ll)>=6 and ll[5] == "-":
            strand_flap = 1
        else:
            strand_flap = 0
        start = int(ll[1])
        end = int(ll[2])
        S = max(0,start - extend)
        E = end + extend
#        S = int(ll[1])
#        E = int(ll[2])
        outdata = ll
        try:
            plus_signal=(plus.summarize(ll[0],S,E,(E-S)))
            minus_signal=(minus.summarize(ll[0],S,E,(E-S)))
            if plus_signal and minus_signal:
                plus_tmp = list(plus_signal.sum_data)
                minus_tmp = list(minus_signal.sum_data)

                if strand_flap == 1:
                    thisdata_tmp = minus_tmp[::-1] + plus_tmp[::-1]#map(round,thisdata_tmp,[4]*(E-S))[::-1]
                else:
                    thisdata_tmp = plus_tmp + minus_tmp
                thisdata = thisdata_tmp#map(round,thisdata_tmp,[4]*len(thisdata_tmp))
        except:
            pass              
        outdata.extend(thisdata)
            # ll.extend(list(signal.sum_data/signal.valid_count))
        outf.write("\t".join(map(str,outdata))+"\n")
    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """calculate the average signal within full len of input peak, multiple bw files separated by comma"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="input bed file")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="output bedEX file")
    optparser.add_option("-p","--plusbw",dest="plusbw",type="str",default = "",
                         help="plus bigwig cuts")
    optparser.add_option("-n","--minusbw",dest="minusbw",type="str",default = "",
                         help="minus bigwig cuts")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",default = "",
                         help="folder of bigwig files")
    optparser.add_option("--ext",dest="extend",type="int",default = 40,
                         help="extend size from center")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    if not inputfile:
        optparser.print_help()
        sys.exit(1)

    get_signal(inputfile,output,options.plusbw,options.minusbw,options.bwfolder,options.extend)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

