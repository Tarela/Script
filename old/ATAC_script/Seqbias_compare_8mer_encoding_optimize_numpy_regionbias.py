#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description: this one seems to have bugs in fetch sequence

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
    pBG = {}
    nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        ll = line.split()
        name = ll[0]
        pBG[name] = float(ll[1])
        nBG[name] = float(ll[2])
    inf.close()
    return pBG,nBG

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
        if seqd[seq] == 0:
            continue
        # constant
        s = [1]
        # single nucleotides 
        for i,elem in enumerate( l ):
            s += code.d1[l[i]] 
        # interaction of adjacent nucleotides 
        for i,elem in enumerate( l[0:-1] ):
            s += code.d2[ l[i]+l[i+1] ]

        mat.append(s)

        y.append(numpy.log(seqd[seq]))
        
    y = numpy.array(y)
    b = ols( numpy.array(mat), y )
    #print b
    b0 = b[0]
    b1 = b[1:(1+3*len(seq))]
    b1 = numpy.reshape(b1,(-1,3))
    b2 = b[(1+3*len(seq)):len(b)]
    b2 = numpy.reshape(b2,(-1,9))
    return b0,b1,b2




def predict2(seq,b0,b1,b2,w):
    """
    prediction based on w-mer using 
    mono- and di-nucleotides
    """
    n = len(seq)
    code = encoding()
    mono = numpy.array( list(seq) )
    di   = numpy.array( [ seq[i]+seq[i+1] for i,elem in enumerate(seq[0:-1]) ] ) 

    # make design matrix from sequence
    # mono-
    X1 = numpy.zeros( (len(mono),3) )
    # di-
    X2 = numpy.zeros( (len(di),9)   )

    for nuc,v in code.d1.iteritems():
        X1[mono==nuc] = numpy.array(v)
    for dinuc,v in code.d2.iteritems(): 
        X2[di==dinuc] = numpy.array(v)

    Y0 = b0*numpy.ones(len(mono))
    Y1 = numpy.dot(X1,numpy.transpose(b1))
    Y2 = numpy.dot(X2,numpy.transpose(b2))

    Y = Y0[0:(n-(w-1))]
    for i in range(w):
        Y += Y1[i:(n-(w-1)+i),i] 

    for i in range(w-1):
        Y += Y2[i:(n-(w-1)+i),i]

    return Y


def sitepro_scan(peak,out_bed,w_plus,w_minus,bg0,tag_span,bias_span,gen,lflank,rflank,offset,bpshift):
 
    t = time.time()
    nmer = lflank + rflank
    genome = twobitreader.TwoBitFile(gen)
    
    pBG,nBG = readBG(bg0)

    code = encoding() 
    b0s0,b1s0,b2s0 = paramest(pBG)

    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    outf = open(out_bed,'w')
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        p_sum = list(w_plus_H.summarize(chrm,start-tag_span,end+tag_span,end-start+2*tag_span).sum_data)
        n_sum = list(w_minus_H.summarize(chrm,start-tag_span,end+tag_span,end-start+2*tag_span).sum_data)

        ### new method
        ### 1. fetch seq
        seq_start_plus = start - tag_span - offset + 1 + bpshift - lflank - bias_span
        seq_start_minus = start - tag_span + 1 - bpshift - rflank - bias_span
        seq_end_plus = end + tag_span - 1 + bpshift + rflank + bias_span - 1
        seq_end_minus = end + tag_span - 1 + offset - bpshift + lflank + bias_span - 1
        
        seq_start = min(seq_start_plus,seq_start_minus)
        if seq_start < 0 :
      #      print line
            continue
        seq_end = max(seq_end_plus,seq_end_minus)
        code_seq = genome[chrm][seq_start:seq_end].upper()
        
        if 'N' in code_seq :
       #     print line 
            continue
        
        ### 2. fetch base
        seq_f = code_seq[ (seq_start_plus - seq_start) : (seq_end_plus - seq_start) ]
        seq_r = code_seq[ (seq_start_minus - seq_start) : (seq_end_minus - seq_start) ]
        
        yf  = predict2(seq_f,b0s0,b1s0,b2s0,8)
        yr  = predict2(revcomp(seq_r),b0s0,b1s0,b2s0,8)[::-1]
        z   = yf + yr
        
        zp = z[(offset-1):]
        zn = z[:(-offset+1)]
        linear_zp = pow(numpy.e,zp)
        linear_zn = pow(numpy.e,zn)
        pnew = list(pow(numpy.e,zp))[bias_span:len(zp)-(bias_span-1)]
        nnew = list(pow(numpy.e,zn))[bias_span:len(zn)-(bias_span-1)]

        pbig_mean = linear_zp[:len(zp)-(2*bias_span-1)]/ (2.0*bias_span)
        pbig_sq = zp[:len(zp)-(2*bias_span-1)]/ (2.0*bias_span)
        nbig_mean = linear_zn[:len(zn)-(2*bias_span-1)]/ (2.0*bias_span)
        nbig_sq = zn[:len(zn)-(2*bias_span-1)]/ (2.0*bias_span)

        for i in range(1,2*bias_span):
            pbig_mean += linear_zp[i:len(zp)-(2*bias_span-1-i)]/ (2.0*bias_span)
            pbig_sq = zp[i:len(zp)-(2*bias_span-1-i)]/ (2.0*bias_span)
            nbig_mean += linear_zn[i:len(zn)-(2*bias_span-1-i)]/ (2.0*bias_span)
            nbig_sq = zn[i:len(zn)-(2*bias_span-1-i)]/ (2.0*bias_span)
        pbig_sq = pow(numpy.e,pbig_sq)
        nbig_sq = pow(numpy.e,nbig_sq)

        ## get predicted seqbias
        pnew_assign =[]
        nnew_assign =[]
        pbig_mean_assign=[]
        nbig_mean_assign=[]
        pbig_sq_assign=[]
        nbig_sq_assign=[]
        pobs = []
        nobs = []

        for bp in range(len(p_sum)- 2*tag_span):
            
            ptotal = 1.0*( sum(p_sum[bp:(bp+tag_span-bias_span)]) + sum(p_sum[(bp+tag_span+bias_span):(bp+2*tag_span)]) ) * 2*tag_span / (2*tag_span - 2*bias_span)
            ntotal = 1.0*( sum(n_sum[bp:(bp+tag_span-bias_span)]) + sum(n_sum[(bp+tag_span+bias_span):(bp+2*tag_span)]) ) * 2*tag_span / (2*tag_span - 2*bias_span)
            
            pobs.append( sum(p_sum[(bp+tag_span-bias_span):(bp+tag_span+bias_span)]) *0.5/bias_span )
            nobs.append( sum(n_sum[(bp+tag_span-bias_span):(bp+tag_span+bias_span)]) *0.5/bias_span )
            
            pnew_assign.append(ptotal * pnew[bp+tag_span]/sum(pnew[bp:(bp+2*tag_span)]))
            nnew_assign.append(ntotal * nnew[bp+tag_span]/sum(nnew[bp:(bp+2*tag_span)]))
                
            pbig_mean_assign.append(ptotal * pbig_mean[bp+tag_span]/sum(pbig_mean[bp:(bp+2*tag_span)]))
            nbig_mean_assign.append(ntotal * nbig_mean[bp+tag_span]/sum(nbig_mean[bp:(bp+2*tag_span)]))
            
            pbig_sq_assign.append(ptotal * pbig_sq[bp+tag_span]/sum(pbig_sq[bp:(bp+2*tag_span)]))
            nbig_sq_assign.append(ntotal * nbig_sq[bp+tag_span]/sum(nbig_sq[bp:(bp+2*tag_span)]))
                
        
        #print len(pobs),len(pnew_assign),len(pbig_mean_assign),len(pbig_sq_assign)                
        ### write  real cleavage , seqbias , seqbias predicted cleavage
        newll = ll + pobs + nobs + pbig_mean_assign + nbig_mean_assign + pbig_sq_assign + nbig_sq_assign# +  p_sum[span:(len(p_sum)-span)] + n_sum[span:(len(n_sum)-span)] +pnew[span:(len(p_sum)-span)] + nnew[span:(len(n_sum)-span)] + pnew_assign + nnew_assign
#        newll = ll + p_sum[span:(len(p_sum)-span)] + n_sum[span:(len(n_sum)-span)] + praw_assign +nraw_assign +px_assign +nx_assign+  pnew_assign + nnew_assign
#        newll = ll + p_sum[span:(len(p_sum)-span)] + n_sum[span:(len(n_sum)-span)] + pnew[span:(len(n_sum)-span)] + nnew[span:(len(n_sum)-span)] + pnew_assign + nnew_assign + pf6_assign + pf8_assign
        outf.write("\t".join(map(str,newll))+"\n")
#    print "predict cut time :",time.time()-t

    outf.close()
    inf.close()
 
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="if raw , -i is 5 column summit200bp file")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/Bam/GM12878_ATACseq_50k_p.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/Bam/GM12878_ATACseq_50k_m.bw",
                         help="")
    optparser.add_option("-b","--bg",dest="bgmatrix",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/biasoffset/up4down4/GM12878_ATACseq_50k_bias_onPeak_u4d4s0_newcode.txt",
                         help="sequence bias matrix at shift 0")
 #   optparser.add_option("--bg2",dest="bgmatrix2",type="str",default = "/mnt/Storage2/home/huse/ATAC/Data/ATAC/biasoffset/up3down3/GM12878_ATACseq_50k_bias_onPeak_u3d3s2_newcode.txt",
 #                        help="sequence bias matrix at shift 2")


#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 50,
                         help="region for get total signal in single bp, default = 50 means +-50bp(total 100bp) signal as total for each bp")
    optparser.add_option("--Bspan",dest="Bspan",type="int",default = 12,
                         help="region for get average bias in single bp, default = 12 means +-12bp(total 24bp) signal as big bias for each bp")

    optparser.add_option("--genome",dest="genome",type="str",default = "/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit",
                         help="2bit format")
    optparser.add_option("--left",dest="leftflank",type="int",default = 4,
                         help="flnaking region for seqbias , 8-mer means left=right=4")
    optparser.add_option("--right",dest="rightflank",type="int",default = 4,
                         help="flnaking region for seqbias , 8-mer means left=right=4")
    optparser.add_option("--offset",dest="offset",type="int",default = 9,
                         help="offset related, distance of pair of +/- related cut,default = 9")
    optparser.add_option("--bpshift",dest="bpshift",type="int",default = 0,
                         help="bp of shift center when calculate bias for given cut site,default is 0 ")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    out_bed = options.output
    w_plus = options.w_plus
    w_minus = options.w_minus
    bgmatrix = options.bgmatrix
    gen = options.genome
    lflank = options.leftflank
    rflank = options.rightflank
    Cspan = options.Cspan
    Bspan = options.Bspan
    offset = options.offset
    sitepro_scan(interval,out_bed,w_plus,w_minus,bgmatrix,Cspan,Bspan,gen,lflank,rflank,offset,options.bpshift)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

