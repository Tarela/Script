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
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------


def make_cut(cutbw_handle,ll,span,flength):
    '''
    flength is fetch length, treat different for +/- strand motif
    return a list
    '''
    total = list(cutbw_handle.summarize(ll[0],max(0,int(ll[1])-flength),int(ll[2])+flength,int(ll[2])-int(ll[1])+ flength *2).sum_data)
    #ts = sum(total)
    if ll[5] != '-':
        out = total[ ( flength - span ) : ( flength + int(ll[2]) - int(ll[1]) + span ) ]
    else:
        out = total[ ( flength - span ) : ( flength + int(ll[2]) - int(ll[1]) + span ) ][::-1]
    return out
    
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
def read_distribution_dict(distribution,pre_uplim,cut_uplim):
    inf = open(distribution)
    dis_dict = {}
    
    for line in inf:
        if line.startswith('predict'):
            continue
        ll = line.split()
        if float(ll[0]) > pre_uplim:
            #print ll[0]
            continue
        raw = numpy.array(map(int,ll[1:(1+cut_uplim)]))
        cum = raw.cumsum()
        dis_dict[float(ll[0])] = cum
    inf.close()
    print 1
    return dis_dict
def maek_sti(cleavage,seqbias):
    c_sum = sum(cleavage)
    s_sum = sum(seqbias)
    stimulate_seqbias = [0]*len(seqbias)
    proportion_seqbias = numpy.array(seqbias)/s_sum
    cumsum_seqbias = proportion_seqbias.cumsum()
    #print seqbias
   # print proportion_seqbias
   # print cumsum_seqbias
    #print c_sum
    for i in range(int(c_sum)):
        rd = random.random()
        judge = rd < cumsum_seqbias
       # print rd
      #  print cumsum_seqbias
      #  print judge
        for j in range(len(judge)) : 
            if judge[j] == True:
                stimulate_seqbias[j] += 1
                break
    return stimulate_seqbias
def getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,gen,left,right,matrix_p,matrix_n,prelim,cutlim,fetch_length=100):
    
    genome = twobitreader.TwoBitFile(gen)
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    dis_plus = (read_distribution_dict(matrix_p,prelim,cutlim))
    dis_minus= (read_distribution_dict(matrix_n,prelim,cutlim))

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
        seq = genome[chrom][(start-pspan-left):(end + pspan+right)]
        pout = make_cut(pcutbw,ll,pspan,fetch_length)
        nout = make_cut(ncutbw,ll,pspan,fetch_length)
        ptotal = sum(pout)
        ntotal = sum(nout)
        if strand == "-":
            pout,nout = nout,pout
        #    Ipout,Inout = Inout,Ipout
        if pout == 'NA':
            continue        

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
        pbgtotal = sum(pbglist)
        nbgtotal = sum(nbglist)
        #pout,nout,pbglist,nbglist
        simulate1_plus = []
        simulate1_minus = []
        simulate2_plus = []
        simulate2_minus = []
        for i in range(len(pout)):
            paraw = (pbglist[i]/pbgtotal)*ptotal
            naraw = (nbglist[i]/nbgtotal)*ntotal
            pa = int(paraw)
            na = int(naraw) 
            if dis_plus.has_key(pa):    
                samplenumberP = random.random()*(dis_plus[pa][-1])
                for realcleavage in range(cutlim):
                    if samplenumberP <= dis_plus[pa][realcleavage]: 
                        si_P = realcleavage
                        break
            else:
                si_P = pa
                print 0
            if dis_minus.has_key(na):
                samplenumberN = random.random()*(dis_minus[na][-1])
                for realcleavage in range(cutlim):
                    if samplenumberN <= dis_minus[na][realcleavage]: 
                        si_N = realcleavage
                        break  
            else:
                si_N = na
                print 0
            #print dis_minus.keys()
            
            simulate1_plus.append(si_P)
            simulate1_minus.append(si_N)
        #print pout,pbglist
        simulate2_plus = maek_sti(pout,pbglist)
        simulate2_minus = maek_sti(nout,nbglist)
                
        newll = ll +  pout + nout + simulate1_plus + simulate1_minus  + simulate2_plus + simulate2_minus + pbglist + nbglist
        outf.write("\t".join(map(str,newll))+"\n")
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

    optparser.add_option("--left",dest="left",type="int",default=3,
                         help="")
    optparser.add_option("--right",dest="right",type="int",default=3,
                         help="")


    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 25,
                         help="default = 25 , means total length = 50")

    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="2bit format")
    optparser.add_option("--m1",dest="m_plus",type="str",
                         help="predict vs observation matrix 1 is plus , 2 is minus")
    optparser.add_option("--m2",dest="m_minus",type="str",
                         help="")
    optparser.add_option("--prelim",dest="prelim",type="int",default = 500,
                         help="preidct cut up limit , any predicted number greater than this will get Zscore as 0")
    optparser.add_option("--cutlim",dest="cutlim",type="int",default = 5000,
                         help="real cut up limit , any predicted number greater than this will get Zscore as 0")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    BGmatrix = options.bgmatrix
    pcut = options.pcut
    ncut = options.ncut
    genome = options.genome
    pspan = options.pspan
    matrix_p = options.m_plus
    matrix_n = options.m_minus
    
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
#    getsignal(inputfile,outputfile,bgmatrix,pcut,ncut,pspan,genome,options.left,options.right,fetch_length=100)
    getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,pspan,genome,options.left,options.right,matrix_p,matrix_n,options.prelim,options.cutlim,fetch_length=100)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
