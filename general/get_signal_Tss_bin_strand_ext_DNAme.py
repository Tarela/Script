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

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
def bwsig_pattern(bwfile,chrm,start,end,points):
    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
    sigPat = sp(cmd)[0].strip().split("\t")
    if len(sigPat) <=1:
        return ['n/a']*points
    else:
        return sigPat

def summary_DNAme_pattern(indata,binsize,strand):
    #print indata
    part = len(indata)/binsize
    NCpG = []
    NunmeCpG = []
    NmeCpG = []
    aveDNAme = []
    for i in range(part):
        partdata = numpy.array(indata[(i*binsize):((i+1)*binsize)])
        filterdata = list( partdata[partdata!='n/a'] ) 
        if filterdata == [''] or len(filterdata) == 0:
            NCpG.append(0)
            NmeCpG.append(0)
            NunmeCpG.append(0)
            aveDNAme.append(-1)
            #return -1 
        else:
            realdata = numpy.array(map(float,filterdata))
            NCpG.append(len(realdata))
            NmeCpG.append(len(realdata[realdata > 0.8]))
            NunmeCpG.append(len(realdata[realdata <= 0.2]))
            aveDNAme.append(numpy.mean(realdata))
    if strand == "+":
        return NCpG + NunmeCpG + NmeCpG + aveDNAme
    else:
        return NCpG[::-1] + NunmeCpG[::-1] + NmeCpG[::-1] + aveDNAme[::-1]

def get_signal(inputfile,output,signalbw,extend,binsize,bwfolder):
    signalbw = signalbw.strip().strip(',').split(',')
    if not bwfolder:
        bwfolder = ""
    if not bwfolder.endswith('/'):
        bwfolder += '/'
    bwHs = []
    for sb in signalbw:
        if "/" in sb:#sb.startswith('/mnt/Storage'):
            bwHs.append(sb)
        else:
            bwHs.append(bwfolder + sb)

    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        if "_" in ll[0]:
            continue
        inputlen = len(ll)
#        print ll
        if ll[5] == "+":
            center = int(ll[1])
        else:
            center = int(ll[2])
        S = center - extend
        E = center + extend
        if S < 0:
            continue
        for bwH in bwHs:#bwHandle:
         #   print ll[0],S,E,N
            
            tmp_result = bwsig_pattern(bwH,ll[0],S,E,E-S)
            result = summary_DNAme_pattern(tmp_result, binsize, ll[5])
            ll.extend(result)
            #print bwH, result
            #if 'n/a' in result or len(result) != N:
            #    continue#print result
            #if ll[5] == "+":
            #    ll.extend(result)
            #else:
            #    ll.extend(result[::-1])
        #if 'n/a' in ll or len(ll) != (inputlen + N*len(signalbw)):
        #    pass
        #else: 
        outf.write("\t".join(map(str,ll))+"\n")
    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")
    optparser.add_option("-w","--signalbw",dest="bw",type="str",default = "",
                         help="")
    optparser.add_option("--ext",dest="ext",type="int",default = "2000",
                         help="")
    optparser.add_option("--binsize",dest="binsize",type="int",default = "200",
                         help="binsize, default is 200bp")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",
                         help="")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    signalbw = options.bw
    bwfolder = options.bwfolder
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,signalbw,options.ext,options.binsize,bwfolder)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

