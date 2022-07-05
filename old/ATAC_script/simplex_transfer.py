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
import numpy
import twobitreader
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def simplex_transfer(biastable,outtable):
    fp = open(biastable)
    
    seqd = {}
    for line in fp.readlines():
        f = line.split()
        try:
            seqd[f[0]] = float(f[1])
        except:
            pass 

    d1  = { 'A': [1,-1,-1], 'C': [-1,1,-1], 'G': [-1,-1,1], 'T':[1,1,1] }
    d2  = { 
    'AA': [+1,-1,-1,-1,+1,+1,-1,+1,+1], 
    'AC': [-1,+1,-1,+1,-1,+1,+1,-1,+1],
    'AG': [-1,-1,+1,+1,+1,-1,+1,+1,-1],
    'AT': [+1,+1,+1,-1,-1,-1,-1,-1,-1],
    'CA': [-1,+1,+1,+1,-1,-1,-1,+1,+1],
    'CC': [+1,-1,+1,-1,+1,-1,+1,-1,+1],
    'CG': [+1,+1,-1,-1,-1,+1,+1,+1,-1], 
    'CT': [-1,-1,-1,+1,+1,+1,-1,-1,-1],
    'GA': [-1,+1,+1,-1,+1,+1,+1,-1,-1],
    'GC': [+1,-1,+1,+1,-1,+1,-1,+1,-1],
    'GG': [+1,+1,-1,+1,+1,-1,-1,-1,+1], 
    'GT': [-1,-1,-1,-1,-1,-1,+1,+1,+1],
    'TA': [+1,-1,-1,+1,-1,-1,+1,-1,-1],
    'TC': [-1,+1,-1,-1,+1,-1,-1,+1,-1],
    'TG': [-1,-1,+1,-1,-1,+1,-1,-1,+1], 
    'TT': [+1,+1,+1,+1,+1,+1,+1,+1,+1] }

    mat,y  = [],[]
    
    ofp = open(outtable,'w')
    
    for j,seq in enumerate(seqd):
        l = list(seq)
        # constant
        s = [1]
        
        # single nucleotides 
        for i,elem in enumerate( l ):
            s += d1[l[i]] 
        # interaction of adjacent nucleotides 
        for i,elem in enumerate( l[0:-1] ):
            s += d2[ l[i]+l[i+1] ]
            
        outll= [seq,numpy.log(seqd[seq])]+s
        outline = "\t".join(map(str,outll))+"\n"
        ofp.write(outline)
    fp.close()
    ofp.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
  
    optparser.add_option("-i","--biastable",dest="biastable",type="str",
                         help="")              
  
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")              
                         
#========minor options=============

    (options,args) = optparser.parse_args()
    biastable = options.biastable
    out = options.out
    if not biastable or not out:
        optparser.print_help()
        sys.exit(1)
    
    simplex_transfer(biastable,out)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)




