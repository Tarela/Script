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
import copy,time

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def read_assign(assign,outp,uplim,c_range):
    p_outlist = []
    n_outlist = []
    for i in range(uplim):
        p_outlist.append([0]*(c_range+1))
        n_outlist.append([0]*(c_range+1))
    inf = open(assign)
    for line in inf:
        ll = line.split()
        each = (len(ll)-5)/6
        p_cut = ll[5:(5+each)]
        n_cut = ll[(5+each):(5+each*2)]
        p_bias = ll[(5+each*2):(5+each*3)]
        n_bias = ll[(5+each*3):(5+each*4)]
        p_assign = ll[(5+each*4):(5+each*5)]
        n_assign = ll[(5+each*5):(5+each*6)]
        for i in range(len(p_assign)):
            pa = float(p_assign[i])
            na = float(n_assign[i])
            pc = float(p_cut[i])
            nc = float(n_cut[i])
            
            p_outilst[int(pa)][int(pc)] += 1
            n_outlist[int(na)][int(pc)] += 1

           #     #print up,pa
           #     if pa < (up + 1) and pa >= up:
           #         #print 1
           #         if pc <= c_range: 
           #             p_outlist[up][int(pc)] += 1
           #             
           #     if na < (up + 1) and na >= up:
           #         if nc <= c_range: 
           #             n_outlist[up][int(nc)] += 1
    inf.close()
    outfp = open(outp,'w')
    #outfn = open(outn,'w')
    title =  ['assign_seqbias']
    for i in range(c_range+1):
        colname = 'assign_'+str(i)
        title.append(colname)
    outfp.write("\t".join(map(str,title))+"\n")
    #outfn.write("\t".join(map(str,title))+"\n")
    
    for i in range(uplim):
        olist = []
        for j in range(len(p_outlist[i])):
            olist.append(p_outlist[i][j]+n_outlist[i][j])
        pl = ["as_"+str(i)+"_"+str(i+1)] + olist
        #nl = ["as_"+str(i)] + n_outlist[i]
        outfp.write("\t".join(map(str,pl))+"\n")
        #outfn.write("\t".join(map(str,nl))+"\n")
    outfp.close()
    outfn.close()
            
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
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")
 #   optparser.add_option("-n","--outn",dest="outn",type="str",
 #                        help="")
                         
    optparser.add_option("--up",dest="uplim",type="int",default=20,
                         help="uplimit of the assigned cleavage ")
    optparser.add_option("--c_range",dest="cleavage_range",type="int",default=500,
                         help="range of x-axis of histogram")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    out = options.out
  #  outn = options.outn
    
    uplim = options.uplim
    c_range = options.cleavage_range
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    read_assign(inputfile,out,uplim,c_range)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


