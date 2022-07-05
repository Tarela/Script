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
import numpy
from numpy import linalg as la
from copy import deepcopy
#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def revcomp(seq):
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A' }
    seqt = list(seq)
    seqt.reverse()
    r = ''.join( [ rc[x] for x in seqt] )
    return r

def ols(X,y):
    # b=(X'X)^-1X'y
    XT = numpy.transpose(X)
    Z = numpy.dot(XT,X) 
    A = numpy.dot( la.inv( numpy.dot(XT,X) ), XT )
    b = numpy.dot(A,y)

#    yp =  numpy.dot(X,b)
    return b

def readBG(bgmatrix):
    #Mgenomepos = 0
    #Nreadscount = 0
    pBG = {}
    inf = open(bgmatrix)
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        name = ll[0]
#        pBG[name] = pow(numpy.e,float(ll[1]))
        pBG[name] = float(ll[2])
    #    Mgenomepos += int(ll[3])
    #    Nreadscount += int(ll[2])        
    inf.close()
    return pBG,len(name)/2#, Nreadscount, Mgenomepos

class encoding(object):
    def __init__(self):
        
        self.mononuc = ['A','C','G','T'] 
        self.dinuc   = ['AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']
 
        self.d1  = { 'A': [1,-1,-1], 'C': [-1,1,-1], 'G': [-1,-1,1], 'T':[1,1,1] }
        self.d2  = { 
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
            #'GC': [+1,-1,+1,+1,-1,+1,+1,-1,+1],
            #'TG': [-1,-1,+1,-1,-1,+1,-1,-1,-1], 

        # parameters
        self.b0 = 0
        self.b1 = []
        self.b2 = []

        self.B1 = numpy.array( [
            [+1,-1,-1],
            [+1,-1,-1],
            [+1,-1,-1],
            [+1,-1,-1],
            [-1,+1,-1], 
            [-1,+1,-1],
            [-1,+1,-1], 
            [-1,+1,-1], 
            [-1,-1,+1], 
            [-1,-1,+1], 
            [-1,-1,+1], 
            [-1,-1,+1], 
            [+1,+1,+1], 
            [+1,+1,+1], 
            [+1,+1,+1],
            [+1,+1,+1]] 
        )

        self.B2 = numpy.array([
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1], 
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1], 
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1], 
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1] ]
        )
 
        self.H = numpy.array( 
            [[+1,-1,-1,-1,+1,+1,-1,+1,+1],
             [-1,+1,-1,+1,-1,+1,+1,-1,+1],
             [-1,-1,+1,+1,+1,-1,+1,+1,-1],
             [+1,+1,+1,-1,-1,-1,-1,-1,-1],
             [-1,+1,+1,+1,-1,-1,-1,+1,+1],
             [+1,-1,+1,-1,+1,-1,+1,-1,+1],
             [+1,+1,-1,-1,-1,+1,+1,+1,-1],
             [-1,-1,-1,+1,+1,+1,-1,-1,-1],
             [-1,+1,+1,-1,+1,+1,+1,-1,-1],
             [+1,-1,+1,+1,-1,+1,-1,+1,-1],
             [+1,+1,-1,+1,+1,-1,-1,-1,+1],
             [-1,-1,-1,-1,-1,-1,+1,+1,+1],
             [+1,-1,-1,+1,-1,-1,+1,-1,-1],
             [-1,+1,-1,-1,+1,-1,-1,+1,-1],
             [-1,-1,+1,-1,-1,+1,-1,-1,+1],
             [+1,+1,+1,+1,+1,+1,+1,+1,+1]] )

    def check(self):
        A = numpy.hstack( (self.B1, numpy.hstack( (self.B2, self.H) ) ))
        #print numpy.dot( numpy.transpose( A ), A ) 

    def rparam(self):
        """
        compute parameters for reverse orientation
        """
        rcmono = {'A':'T','C':'G','G':'C','T':'A'} 
        rcdi   = {'AA':'TT', 'AC':'GT', 'AG':'CT', 'AT':'AT', 'CA':'TG', 'CC':'GG', 'CG':'CG', 'CT':'AG', 'GA':'TC', 'GC':'GC', 'GG':'CC', 'GT':'AC', 'TA':'TA', 'TC':'GA', 'TG':'CA', 'TT':'AA'}
         
        A = numpy.hstack( ( numpy.transpose(numpy.array([16*[1]])), numpy.hstack( (self.B1, numpy.hstack( (self.B2, self.H) ) ))) )
        #print numpy.dot( numpy.transpose( A ), A ) 

        R = []
        for i,elem in enumerate(self.dinuc):
            rc = rcdi[elem]
            s = [1] + self.d1[rc[0]] + self.d1[rc[1]] + self.d2[rc]   
            R.append(s)  
        #print  numpy.array(R) 
        #print numpy.dot( numpy.transpose( R ), R ) 
        self.RCMAP =  numpy.dot( 1.0/16*numpy.transpose( A ), R ) 


def paramest(seqd):
    """
    Implementation of
    Maximally Efficient Modeling of DNA Sequence Motifs at All Levels of Complexity
    Gary D. Stormo
    Genetics, 2011
    """

    code=encoding()
    #code.check()
    #code.rparam()
    mat,y  = [],[]
    for j,seq in enumerate(seqd):
        l = list(seq)    
        #if seqd[seq] <= 0:
        #    continue
        # constant
        s = [1]
        # single nucleotides 
        for i,elem in enumerate( l ):
            s += code.d1[l[i]] 
        # interaction of adjacent nucleotides 
        for i,elem in enumerate( l[0:-1] ):
            s += code.d2[ l[i]+l[i+1] ]

        mat.append(s)

        y.append((seqd[seq]))
#        y.append(numpy.log(seqd[seq]))
#        if seq == "ACTCGCAA":
#            print numpy.log(seqd[seq])
        
    y = numpy.array(y)
    #print list(y)
    b = ols( numpy.array(mat), y )
    #print b
    b0 = b[0]
    b1 = b[1:(1+3*len(seq))]
    #b1 = numpy.reshape(b1,(-1,3))
    b2 = b[(1+3*len(seq)):len(b)]
    #b2 = numpy.reshape(b2,(-1,9))
    return b,b0,b1,b2


def seq2biasParm(useseq,code):
    #param_length_b0 = 1
    #param_length_b1 = len(seq)*3 
    #param_length_b2 = (len(seq)-1)*9

#    item_b0 = 0
#    item_b1 = numpy.array([0]*param_length_b1)
#    item_b2 = numpy.array([0]*param_length_b2)

    item_b1 = []
    item_b2 = []

    for n in range(len(useseq)):
        bp1 = useseq[n]
        simplex = code.d1[bp1]
        item_b1.extend(simplex)
    for n2 in range(len(useseq) -1):
        bp2 = useseq[n2:(n2+2)]
        simplex = code.d2[bp2]
        item_b2.extend(simplex)

    seq_item = numpy.array([1] + item_b1 + item_b2)

    return seq_item 




def get_regionLevel_simplex_parameters(inputbed,outputbed,flank,genome2bit,ext):
    simplex_code = encoding()
    #biasdict,flank = readBG(biasmat)
    #B,B0,B1,B2 = paramest(biasdict)
    flank = int(flank)
#    outitem = seq2biasParm("ACTCGCAA",B,simplex_code)
    #print B
    genome = twobitreader.TwoBitFile(genome2bit)
#    seq = genome[chrm][(int(ll[1])-flank):(int(ll[1])+flank)].upper()

    inf = open(inputbed)
    outf = open(outputbed,'w')
    simplex_len = len(list(seq2biasParm("A"*(flank*2),simplex_code)))
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        if "_" in chrm:
            continue
        center = (int(ll[1]) + int(ll[2]))/2
        start = max(0,center-ext)
        end = center + ext

        plusSequence = genome[chrm][(start-flank):(end+flank)].upper()
        minusSequence = genome[chrm][(start-flank+1):(end+flank+1)].upper()

        #plus_data = "NA"
        #minus_data = "NA"
        plus_data = numpy.array([0.0]*simplex_len)
        minus_data = numpy.array([0.0]*simplex_len)
        for i in range(len(plusSequence)-2*flank):
            pseq = plusSequence[i:(i+2*flank)].upper()
            if not "N" in pseq:
                p_out = seq2biasParm(pseq,simplex_code)
        #        if plus_data == "NA":
        #            plus_data = numpy.array([0.0]*len(p_out))
        #            minus_data = numpy.array([0.0]*len(p_out))
                plus_data += p_out

        for i in range(len(minusSequence)-2*flank):
            tmpseq = minusSequence[i:(i+2*flank)]
            if not "N" in tmpseq:
                mseq = revcomp(tmpseq).upper()
                m_out = seq2biasParm(mseq,simplex_code)
                minus_data += m_out

        newll = ll + list(plus_data) + list(minus_data)
        outf.write("\t".join(map(str,newll))+"\n")

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
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="inputbed file")
    optparser.add_option("-o","--outputbed",dest="outputbed",type="str",
                         help="output bed file, simplex item for both +strand, - strand, appended in the tail of inputbed")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/scratch/sh8tv/Data/Genome/hg38/hg38.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=4,
                         help="flanking region for n-mer , default =4 means n=8")
    optparser.add_option("-e","--ext",dest="ext",type="int",default=100,
                         help="extend length from center")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputbed = options.outputbed
    flank = options.flank

    if not inputbed or not outputbed:
        optparser.print_help()
        sys.exit(1)
    get_regionLevel_simplex_parameters(inputbed,outputbed,flank,options.sequence,options.ext)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


