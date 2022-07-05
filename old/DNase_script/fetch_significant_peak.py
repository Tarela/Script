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
import string
from CistromeAP.taolib.CoreLib.BasicStat.Func import *
from CistromeAP.jianlib.BwReader import BwIO
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit() 
import math
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def get_signal(inputfile,out1,out2,span):

    inf = open(inputfile)
    outf1 = open(out1,'w')
    outf2 = open(out2,'w')
    for line in inf:
        ll = line.split()
        each = 400
        peak = ll[:5]
        #p_cut = ll[5:(5+each)]
        #n_cut = ll[(5+each):(5+each*2)]
        #p_bias = ll[(5+each*2):(5+each*3)]
        #n_bias = ll[(5+each*3):(5+each*4)]
        #p_assign = ll[(5+each*4):(5+each*5)]
        #n_assign = ll[(5+each*5):(5+each*6)]
        #p_zscore = map(float,ll[(5+each*6):(5+each*7)])
        #n_zscore = map(float,ll[(5+each*7):(5+each*8)])
        Pval = ll[(5+each*10 ): (5+each*10 + 388)]
        best_pval = 1
        
        newPval = [1]*span + Pval + [1]*span
        
        on = 0
        if len(newPval) != 400:
           # print len(ll),len(newPval),newPval
            print ll
        for i in range(400):
            nowP = float(newPval[i])
            
            if nowP < 0.01:
                if on == 0:
                    on = 1
                    start = int(ll[1]) + i
                else:
                    pass
            else:
                if on == 0 :
                    pass
                else:
                    on = 0
                    end = int(ll[1])+i
                    outf1.write("\t".join(map(str,[ll[0],start,end]))+"\n")
    
            
            if nowP < best_pval:
                best_pval = nowP
                best_pos = [ll[0],int(ll[1])+i,int(ll[1])+i+1]
        if on == 1:
            end = int(ll[1])+i
            outf1.write("\t".join(map(str,[ll[0],start,end]))+"\n")

        if best_pval < 0.01:
            outf2.write("\t".join(map(str,best_pos))+"\n")

    inf.close()
    outf1.close()
    outf2.close()



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
    optparser.add_option("-o","--out1",dest="out1",type="str",default = "",
                         help="")
    optparser.add_option("-p","--out2",dest="out2",type="str",default = "",
                         help="")
    optparser.add_option("-s","--span",dest="span",type="int",default = 25,
                         help="")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    out2 = options.out2
    out1 = options.out1
    span = options.span
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,out1,out2,span)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


