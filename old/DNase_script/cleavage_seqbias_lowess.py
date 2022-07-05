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
import math
from cistrome import regions as c
import twobitreader
import numpy
import random
#from CistromeAP.taolib.CoreLib.BasicStat.Func import *
#from CistromeAP.jianlib.BwReader import BwIO
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------



    
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
def getZ(p_ob,pa,pred_size,pred_prob,coeff_size,coeff_logprob,nbinom,sdnorm):
    if pa == -1 :
        return 10
        #print p_ob,pa
        #sys.exit(1)
    if  pred_size.has_key(pa)  :
         pSIZE = pred_size[pa]
         pPROB = pred_prob[pa]
    else:
        if pa < 21:
            print pa
        pSIZE = coeff_size[0] + coeff_size[1]*pa
        pPROB = pow(math.e, coeff_logprob[0] + coeff_logprob[1]*math.log(pa) )
    
    if p_ob == 0:
        pval = 1-(nbinom.cdf(p_ob,pSIZE,pPROB)/2)
    else:
        pval = 1 - ( nbinom.cdf(p_ob,pSIZE,pPROB)/2 + nbinom.cdf(p_ob-1,pSIZE,pPROB)/2  )
   # pval = 1-(nbinom.cdf(p_ob,pSIZE,pPROB) - nbinom.cdf(0,pSIZE,pPROB)/2)
    if pval == 1:
        pZscore = -8
    elif pval == 0:
        pZscore = 8
    else:
        pZscore =  sdnorm.isf(pval) 
    if pZscore > 8 :
        print pZscore
    #print p_ob,pa,pval,pZscore
    return pZscore         

def nbinom_zscore(cut,n,p,nbinom,sdnorm):
    pval = 1-nbinom.cdf(cut,n,p)
    Zscore = sdnorm.isf(pval) 
    return Zscore
    

def read_param(parameter):
    pred_size = {}
    pred_prob = {}
    inf = open(parameter)
    for line  in inf:
        if line.startswith('predict'):
            continue
        ll = line.split()
        pred_size[float(ll[0])] = float(ll[1])
        pred_prob[float(ll[0])] = float(ll[2])
    inf.close()
    return pred_size,pred_prob
def read_coeff(coeff):
    
    inf = open(coeff)
    for line in inf:
        ll = line.split()
        if ll[0] == "SIZE":
            coeff_size = [float(ll[1]),float(ll[2])]
        if ll[0] == "logPROB":
            coeff_logprob = [float(ll[1]),float(ll[2])]
    inf.close()
    return coeff_size,coeff_logprob
    
def read_lowess_table(lowess,bgdict):
    inf  = open(lowess)
    lowess_dict = {}
    for i in bgdict.keys():
        bias = bgdict[i]
        lowess_dict[bgdict[i]] = {}
    for line in inf:
        if line.startswith('minbias'):
            continue
        ll = line.split() ## minbias , maxbias , perbp_count , lowess_median_obs
        for bias in lowess_dict.keys():
            if bias >= float(ll[0]) and bias <= float(ll[1]):
                lowess_dict[bias][float(ll[2])] = float(ll[3])
    inf.close()
    return lowess_dict


def getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,lowess,pspan,tspan,gen,left,right):
    pBG,nBG = readBG(BGmatrix)
    lowess_dict = read_lowess_table(lowess,pBG)
    print len(lowess_dict.keys())
    #print sorted(lowess_dict.keys())
    #print lowess_dict[3.62053200333]
    #sys.exit(1)
    # tspan = region get tag
    genome = twobitreader.TwoBitFile(gen)
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))

    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    pspan = pspan - ml/2
    inf.seek(0)
    pBG,nBG = readBG(BGmatrix)
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()

        chrom = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        strand = ll[5]
        
        if int(ll[1])-pspan-tspan < 0 :
            continue
        pout = list(pcutbw.summarize(ll[0],max(0,int(ll[1])-pspan-tspan),int(ll[2])+ pspan + tspan,int(ll[2])-int(ll[1])+ pspan *2 + tspan*2).sum_data)
        nout = list(ncutbw.summarize(ll[0],max(0,int(ll[1])-pspan-tspan),int(ll[2])+ pspan + tspan,int(ll[2])-int(ll[1])+ pspan *2 + tspan*2).sum_data)
        if strand != '-':
            pass
        else:
            pout,nout = nout[::-1],pout[::-1]
        seq = genome[chrom][(start-pspan-tspan-left):(end + pspan + tspan+right)]
        ### get seqbias
        if 'N' in seq.upper():
            continue
        #print 1
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - left-right):
            p.append(pBG[pseq[k:k+left+right].upper()])
            n.append(nBG[nseq[k:k+left+right].upper()])
        if strand != '-':
            pbglist = p
            nbglist = n
        else:
            pbglist = n[::-1]
            nbglist = p[::-1]
        
        if len(pbglist) != len(pout):
            print 'len pbglist != len pout , should have bugs'
            sys.exit(1)
            
        p_cleavage = []
        n_cleavage = []
        p_bias =   []
        n_bias = []
        p_lowess = []
        n_lowess = []
 
        for i in range(len(pout)-2*tspan):
            p_ob = pout[i+tspan]
            n_ob = nout[i+tspan]
            ptotal = sum(pout[i:(i+tspan*2)])
            ntotal = sum(nout[i:(i+tspan*2)])
            pbias = pbglist[i+tspan]
            nbias = nbglist[i+tspan]
            ## make lowess 

            p_count_table = lowess_dict[pbias]
            n_count_table = lowess_dict[nbias]
            #print pbias,lowess_dict[pbias].keys()
            #print nbias,lowess_dict[nbias].keys()
            
            pcount = min(int(float(ptotal)/10)/5.0,max(lowess_dict[pbias].keys()))
            ncount = min(int(float(ntotal)/10)/5.0,max(lowess_dict[nbias].keys()))
            plowess = p_count_table[pcount]
            nlowess = n_count_table[ncount]
            
            p_cleavage.append(p_ob)
            n_cleavage.append(n_ob)
            p_bias.append(pbias)
            n_bias.append(nbias)
            p_lowess.append(plowess)
            n_lowess.append(nlowess)

        newll = ll + p_cleavage +  n_cleavage + p_bias +  n_bias+ p_lowess + n_lowess
        outf.write("\t".join(map(str,newll)) + "\n")
    inf.close()
    outf.close()
        
        #print pa
       #print len(pbglist),len(nbglist), len(pout),len(nout)
        


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
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")                         
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="")
    optparser.add_option("-p","--positive_cut",dest="pcut",type="str",
                         help="")
    optparser.add_option("-n","--negative_cut",dest="ncut",type="str",
                         help="")
    optparser.add_option("-l","--lowess_table",dest="lowess",type="str",
                         help="")
                         
    optparser.add_option("--left",dest="left",type="int",default=3,
                         help="")
    optparser.add_option("--right",dest="right",type="int",default=3,
                         help="")


    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 25,
                         help="default = 25 , means total plot span length = 50")
    optparser.add_option("-t","--tag_span",dest="tspan",type="int",default = 25,
                         help="default = 25 , means total tag span length = 50")

    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="2bit format")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    BGmatrix = options.bgmatrix
    pcut = options.pcut
    ncut = options.ncut
    lowess = options.lowess
    genome = options.genome
    pspan = options.pspan
    tspan = options.tspan

    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
#    getsignal(inputfile,outputfile,bgmatrix,pcut,ncut,pspan,genome,options.left,options.right,fetch_length=100)
    #getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,genome,options.left,options.right,matrix_p,matrix_n,options.prelim,options.cutlim,fetch_length=100)
    getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,lowess,pspan,tspan,genome,options.left,options.right)
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
