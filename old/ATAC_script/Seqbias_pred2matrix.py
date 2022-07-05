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
    draw = make_dict(assign_limit,cut_limit)
    dx = make_dict(assign_limit,cut_limit)
    dsqx = make_dict(assign_limit,cut_limit)
    dmax = make_dict(assign_limit,cut_limit)
    outall = open(outname+"_summary.txt",'w')
    for line in inf:
        ll = line.split()

        pc =   ll[(400*0 + 5):(400*1+5)]
        nc =   ll[(400*1 + 5):(400*2+5)]
        pb =   ll[(400*2 + 5):(400*3+5)]
        nb =   ll[(400*3 + 5):(400*4+5)]
        praw = ll[(400*4 + 5):(400*5+5)]
        nraw = ll[(400*5 + 5):(400*6+5)]
        px  =  ll[(400*6 + 5):(400*7+5)]
        nx  =  ll[(400*7 + 5):(400*8+5)]
        psqx = ll[(400*8 + 5):(400*9+5)]
        nsqx = ll[(400*9 + 5):(400*10+5)]
        pmax = ll[(400*10 + 5):(400*11+5)]
        nmax = ll[(400*11 + 5):(400*12+5)]
        for i in range(400):
            pcut = int(float(pc[i]))
            newll = [pcut,praw[i],px[i],psqx[i],pmax[i]]
            outall.write("\t".join(map(str,newll))+"\n")
            if pcut > cut_limit :
                continue
            if draw.has_key(conver_assign(praw[i])):
                draw[conver_assign(praw[i])][pcut] += 1
            if dx.has_key(conver_assign(px[i])):
                dx[conver_assign(px[i])][pcut] += 1
            if dsqx.has_key(conver_assign(psqx[i])):
                dsqx[conver_assign(psqx[i])][pcut] += 1
            if dmax.has_key(conver_assign(pmax[i])):
                dmax[conver_assign(pmax[i])][pcut] += 1
        for i in range(400):
            ncut = int(float(nc[i]))
            newll = [ncut,nraw[i],nx[i],nsqx[i],nmax[i]]
            outall.write("\t".join(map(str,newll))+"\n")
            if ncut > cut_limit :
                continue
            if draw.has_key(conver_assign(nraw[i])):
                draw[conver_assign(nraw[i])][ncut] += 1
            if dx.has_key(conver_assign(nx[i])):
                dx[conver_assign(nx[i])][ncut] += 1
            if dsqx.has_key(conver_assign(nsqx[i])):
                dsqx[conver_assign(nsqx[i])][ncut] += 1
            if dmax.has_key(conver_assign(nmax[i])):
                dmax[conver_assign(nmax[i])][ncut] += 1
    outraw = open(outname+"_raw.txt",'w')
    outx = open(outname+"_x.txt",'w')
    outsqx = open(outname+"_sqx.txt",'w')
    outmax = open(outname+"_max.txt",'w')
    for assign in sorted(draw.keys()):

        outraw.write("\t".join(map(str,[assign] + draw[assign]))+"\n")
        outx.write("\t".join(map(str,[assign]   + dx[assign]))+"\n")
        outsqx.write("\t".join(map(str,[assign] + dsqx[assign]))+"\n")
        outmax.write("\t".join(map(str,[assign] + dmax[assign]))+"\n")
    outraw.close()
    outx.close()
    outsqx.close()
    outmax.close()
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
                         help="output name, will output 4 table, praw - pc ; px - pc ; psqx - pc ; pmax - pc;")
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

