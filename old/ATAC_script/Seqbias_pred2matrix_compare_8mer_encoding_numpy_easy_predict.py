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
def pred2matrix(inputfile,outname,assign_limit,cut_limit):
    inf = open(inputfile)
#    draw = make_dict(assign_limit,cut_limit)
#    dx = make_dict(assign_limit,cut_limit)
#    dnew = make_dict(assign_limit,cut_limit)
#    dnewlin = make_dict(assign_limit,cut_limit)
    cut_raw_addup = 0
    cut_x_addup = 0
    cut_new_addup = 0    
    peaknum = 0
    for line in inf:
        peaknum += 1
        ll = line.split()

        pc =   ll[(400*0 + 5):(400*1+5)]
        nc =   ll[(400*1 + 5):(400*2+5)]
        px =   ll[(400*2 + 5):(400*3+5)]
        nx =   ll[(400*3 + 5):(400*4+5)]
        pnew = ll[(400*4 + 5):(400*5+5)]
        nnew = ll[(400*5 + 5):(400*6+5)]
        #pnew  =  ll[(400*6 + 5):(400*7+5)]
        #nnew  =  ll[(400*7 + 5):(400*8+5)]
 #       pbnew = ll[(400*8 + 5):(400*9+5)]
 #       nbnew = ll[(400*9 + 5):(400*10+5)]
 #       pnew = ll[(400*10 + 5):(400*11+5)]
 #       nnew = ll[(400*11 + 5):(400*12+5)]
 #       pnewlin = ll[(400*2 + 5):(400*3+5)]
 #       nnewlin = ll[(400*3 + 5):(400*4+5)]
        
        for i in range(400):
       #     cut_raw_addup += pow( int(float(pc[i])) - float(praw[i]) ,2)
            cut_x_addup += pow( int(float(pc[i])) - float(px[i]) ,2)
            cut_new_addup += pow( int(float(pc[i])) - float(pnew[i]) ,2)
        for i in range(400):
       #     cut_raw_addup += pow( int(float(nc[i])) - float(nraw[i]) ,2)        
            cut_x_addup += pow( int(float(nc[i])) - float(nx[i]) ,2)
            cut_new_addup += pow( int(float(nc[i])) - float(nnew[i]) ,2)
           # outf2.write("\t".join([nc[i],nraw[i],nx[i],nnew[i]])+"\n")

    #cut_raw_rmse = pow(cut_raw_addup/(peaknum*800),0.5)
    cut_x_rmse = pow(cut_x_addup/(peaknum*800),0.5)
    cut_new_rmse = pow(cut_new_addup/(peaknum*800),0.5)
    outall = open(outname,'w')
    #outall.write("\t".join(["cut_raw_rmse",str(cut_raw_rmse)])+"\n")
    outall.write("\t".join(["cut_x_rmse",str(cut_x_rmse)])+"\n")
    outall.write("\t".join(["cut_new_rmse",str(cut_new_rmse)])+"\n")
    outall.close()
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
    optparser.add_option("--assignlimit",dest="assignlimit",type="int",default=100,
                         help="range of predicted cut")
    optparser.add_option("--cutlimit",dest="cutlimit",type="int",default=200,
                         help="range of real cut")

#========minor options=============


    (options,args) = optparser.parse_args()

    if not options.inputfile:
        optparser.print_help()
        sys.exit(1)

    inputfile = options.inputfile
    outname = options.outname
    pred2matrix(inputfile,outname,options.assignlimit,options.cutlimit)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

