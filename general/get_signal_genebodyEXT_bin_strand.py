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
        bwfolder = "./"
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
        #center = (int(ll[1]) + int(ll[2]))/2
        #S = max(0,center - extend)
        #E = center + extend
        #C = (int(ll[1]) + int(ll[2]) ) /2
        #S = C - extend
        #E = C + extend
        S=int(ll[1])
        E=int(ll[2])

        for bwHandle in bwHs:
            try:
                signal1=(bwHandle.summarize(ll[0],max(0,S-extend),S,20))
                signal2=(bwHandle.summarize(ll[0],S,E,20))
                signal3=(bwHandle.summarize(ll[0],E,E+extend,20))
                binlen1 = extend*1.0/20
                binlen2 = (E-S)*1.0/20
                binlen3 = extend*1.0/20
                if type(signal1.sum_data) == None or type(signal2.sum_data) == None or type(signal3.sum_data) == None:
                    addsig = [0]*60
                else:
                    addsig1 = signal1.sum_data/binlen1 #float(signal.sum_data/signal.valid_count)
                    addsig2 = signal2.sum_data/binlen2
                    addsig3 = signal3.sum_data/binlen3
                    addsig = list(addsig1) + list(addsig2) + list(addsig3)
            except:
                #print 'c2',line
                addsig = [0]*60#'nan'
               # ll.extend(list(signal.sum_data/signal.valid_count))
            if len(ll)>=6 and ll[5] == "-":
                ll.extend(addsig[::-1])
            else:
                ll.extend(addsig)
        outf.write("\t".join(map(str,ll))+"\n")
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
    optparser.add_option("--ext",dest="ext",type="int",default=10000,
                         help="ext size from center,default=10000")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    bws = options.bws
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,bws,options.bwfolder,options.ext)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


