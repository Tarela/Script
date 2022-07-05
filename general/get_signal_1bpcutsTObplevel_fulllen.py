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

def get_signal(inputfile,output,bwfiles,bwfolder,extend):
    signalbw = bwfiles.strip().strip(',').split(',')

    if not bwfolder:
        bwfolder = ""
    if not bwfolder.endswith('/'):
        bwfolder += '/'

    bwHs = []
    for sb in signalbw:
        if sb.startswith('/'):
            bwHs.append(BigWigFile(open(sb, 'rb')))
        else:
            bwHs.append(BigWigFile(open(bwfolder + sb, 'rb')))

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
#        center = (int(ll[1]) + int(ll[2]))/2
#        S = max(0,center - extend)
#        E = center + extend
        S = int(ll[1])
        E = int(ll[2])
        outdata = []
        for bwHandle in bwHs:
            try:
                signal=(bwHandle.summarize(ll[0],S,E,(E-S)))
                if signal:
                    thisdata_tmp = list(signal.sum_data)
#                    if strand_flap == 1:
#                        thisdata = map(round,thisdata_tmp,[4]*(E-S))[::-1]
#                    else:
                    thisdata = map(round,thisdata_tmp,[4]*(E-S))
                else:
                    thisdata = ["NA"]*(E-S)
            except:
                thisdata = ["NA"]*(E-S)               
            outdata.append(thisdata)
                # ll.extend(list(signal.sum_data/signal.valid_count))

        for pos in range(len(outdata[0])):
            newll = [ll[0],S+pos,S+pos+1]
            for dataorder in range(len(outdata)):
                newll.append(outdata[dataorder][pos])
            outf.write("\t".join(map(str,newll))+"\n")
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
    optparser.add_option("-w","--bws",dest="bws",type="str",default = "",
                         help="bigwig files, comma separate")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",default = "",
                         help="folder of bigwig files")
    optparser.add_option("--ext",dest="extend",type="int",default = 200,
                         help="folder of bigwig files")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    bws = options.bws
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,bws,options.bwfolder,options.extend)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

