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
import math,time
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader

#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def conver_assign(paraw):
    paraw = float(paraw)
    if paraw == 0:
        pa = -1
    elif paraw < 1.0/32:
        pa = 0.0
    elif paraw < 1.0/16:
        pa = 1.0/32
    elif paraw < 1.0/8:
        pa = 1.0/16
    elif paraw < 1.0/4:
        pa = 1.0/8
    elif paraw < 1.0/2:
        pa = 1.0/4
    elif paraw < 1.0:
        pa = 1.0/2
    else:
        pa = int(paraw)
    return pa

def make_dict(assign_limit,cut_limit):
    new_dict = {}
    new_dict[-1] = [0]*(cut_limit+1)
    new_dict[1.0/32] = [0]*(cut_limit+1)
    new_dict[1.0/16] = [0]*(cut_limit+1)
    new_dict[1.0/8] = [0]*(cut_limit+1)
    new_dict[1.0/4] = [0]*(cut_limit+1)
    new_dict[1.0/2] = [0]*(cut_limit+1)
    
    for i in range(assign_limit+1):
        new_dict[i] = [0]*(cut_limit+1)
    return new_dict
def pred2matrix(inputfile,outname,chromlen,PN,start):
    chrdict = {}
    inf = open(chromlen)
    for line in inf:
        ll = line.split()
        chrdict[ll[0]] = int(ll[1])
    inf.close()

    inf = open(inputfile)
    outf = open(outname,'w')
    
    last_chrm = 'NA'
    last_start = 'NA'
    last_end = 'NA'
    last_value = 'NA'

    for line in inf:
        ll = line.split()
        peak_chrm = ll[0]
        peak_start = int(ll[1])
        peak_end = int(ll[2])
        if PN == 'p':
            pnew  =  ll[(400*start + 5):(400*(start+1)+5)]
        else:
            pnew = ll[(400*(start+1) + 5):(400*(start+2)+5)]
#        nnew  =  ll[(400*7 + 5):(400*8+5)]

        if last_chrm == 'NA' :#new file start
        #    newll = [peak_chrm,0,peak_start,-1]
        #    outf.write("\t".join(map(str,newll))+"\n")

            last_chrm = peak_chrm
            last_start = 0
            last_end = peak_start
            last_value = -1
            bpstart = 0### start position from every 400bp
        elif last_chrm != peak_chrm :  ### new chromsome start
            newll = [last_chrm,last_start,last_end,last_value]
            outf.write("\t".join(map(str,newll))+"\n")
            newll = [last_chrm,last_end,chrdict[last_chrm],-1]
            outf.write("\t".join(map(str,newll))+"\n")
            #newll = [this_chrm,0,this_start,-1]
            #outf.write("\t".join(map(str,newll))+"\n")

            last_chrm = peak_chrm
            last_start = 0
            last_end = peak_start
            last_value = -1 
            bpstart = 0  #         
        else: ## same chromsome
            if peak_start <= last_end:## 2 peak overlap
                bpstart = last_end - peak_start
            else:
                newll = [last_chrm,last_start,last_end,last_value]
                outf.write("\t".join(map(str,newll))+"\n")
                last_start = last_end
                last_end = peak_start
                last_value = -1
                bpstart = 0
        
        for pos in range(bpstart,400):
            this_value = pnew[pos]
            
            if this_value != last_value:
                newll = [last_chrm,last_start,last_end,last_value]
                outf.write("\t".join(map(str,newll))+"\n")
                last_value = this_value
                last_start = last_end
                last_end  = last_start + 1
            else:
                last_end += 1
    newll = [last_chrm,last_start,last_end,last_value]
    outf.write("\t".join(map(str,newll))+"\n")
    newll = [last_chrm,last_end,chrdict[last_chrm],-1]
    outf.write("\t".join(map(str,newll))+"\n")
                    
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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="Seqbias_prediction.py result")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="output name, will output 4 table, praw - pc ; px - pc ; pnew - pc ; pnewlin - pc;")
    optparser.add_option("-c","--chromlen",dest="chromlen",type="str",default="/mnt/Storage/data/sync_cistrome_lib/chromLen/hg19.len",
                         help="")
    optparser.add_option("--start",dest="start",type="int",default=6,
                         help="which column to be conform to bigwig, default = 6 means 4rd block,0,2,4,6")

#========minor options=============


    (options,args) = optparser.parse_args()

    if not options.inputfile:
        optparser.print_help()
        sys.exit(1)

    inputfile = options.inputfile
    outname = options.outname
    pred2matrix(inputfile,outname+"_plus.bdg",options.chromlen,'p',options.start)
    pred2matrix(inputfile,outname+"_minus.bdg",options.chromlen,'n',options.start)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

