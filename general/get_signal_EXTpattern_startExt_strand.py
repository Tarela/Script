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

def get_signal(inputfile,output,bwfiles,extend,N,bwfolder):
    signalbw = bwfiles.strip().strip(',').split(',')

    if not bwfolder:
        bwfolder = "./"
    if not bwfolder.endswith('/') and not bwfolder != "":
        bwfolder += '/'

    bwHs = []
    for sb in signalbw:
        if sb.startswith('/') or startswith("./") or startswith("../"):
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
            start = int(ll[2])
            strand_flap = 1
        else:
            start = int(ll[1])
            strand_flap = 0
        S = max(0,start - extend)
        E = start + extend
#        S = int(ll[1])
#        E = int(ll[2])
        outdata = ll
        for bwHandle in bwHs:
            try:
                signal=(bwHandle.summarize(ll[0],S,E,N))
                binlen = (E-S)*1.0/N
                if type(signal.sum_data) == None:
                    print 'c1',line
                    addsig = ["na"]*N
                else:

                    addsig = list(signal.sum_data*1.0/(binlen))#float(signal.sum_data/signal.valid_count)
            except:
                print 'c2',line
                addsig = ["na"]*N#'nan'
               # ll.extend(list(signal.sum_data/signal.valid_count))
            if strand_flap == 1:
                ll.extend(addsig[::-1])
            else:
                ll.extend(addsig)

            # ll.extend(list(signal.sum_data/signal.valid_count))
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
    optparser.add_option("-w","--signalbw",dest="bw",type="str",default = "",
                         help="")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",default = "",
                         help="folder of bigwig files")
    optparser.add_option("--ext",dest="extend",type="int",default = 500,
                         help="extend size from center")
    optparser.add_option("-n","--binnum",dest="N",type="int",default = 100,
                         help="extend size from center")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    if not inputfile:
        optparser.print_help()
        sys.exit(1)

    get_signal(inputfile,output,options.bw,options.extend,options.N,options.bwfolder)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

