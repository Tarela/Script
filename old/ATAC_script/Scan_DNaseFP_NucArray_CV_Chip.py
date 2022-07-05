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
def makeTCFOS(cutbw_plus,cutbw_minus,ll,tspan,ml):
    
    pcut = list(cutbw_plus.summarize(ll[0], ( int(ll[1]) + int(ll[2]) )/2 -tspan ,( int(ll[1]) + int(ll[2]) )/2 +tspan ,2*tspan).sum_data)
    ncut = list(cutbw_minus.summarize(ll[0], ( int(ll[1]) + int(ll[2]) )/2 -tspan ,( int(ll[1]) + int(ll[2]) )/2 +tspan ,2*tspan).sum_data)
    TC =  sum(pcut)+sum(ncut)

    C = sum(pcut[(tspan-ml/2) : (tspan-ml/2+ml)]) + sum(ncut[(tspan-ml/2) : (tspan-ml/2+ml)])
    L = sum(pcut[(tspan-ml/2-ml):(tspan-ml/2)]) + sum(ncut[(tspan-ml/2-ml):(tspan-ml/2)])
    R = sum(pcut[(tspan-ml/2+ml):(tspan-ml/2+2*ml)]) + sum(ncut[(tspan-ml/2+ml):(tspan-ml/2+2*ml)])
 
    FOS = ( (C+1)/(R+1) + (C+1)/(L+1) )
    return TC,FOS

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
    
def make_chip(chipbw_handle,ll,number):
    center = (int(ll[1])+int(ll[2]))/2
    total = (chipbw_handle.summarize(ll[0],max(0,center-100),center+100,number).sum_data)
    vd = (chipbw_handle.summarize(ll[0],max(0,center-100),center+100,number).valid_count)
    outraw = total/(vd+0.0001)
    if ll[5] != '-':
        out = list(outraw)
    else:
        out = list(outraw)[::-1]
    return out
    
def make_mnase(mnasebw_handle,ll,span):
    number = 200
    center = (int(ll[1])+int(ll[2]))/2
    total = (mnasebw_handle.summarize(ll[0],max(0,center-span),center+span,span*2/10).sum_data)
    vd = (mnasebw_handle.summarize(ll[0],max(0,center-span),center+span,span*2/10).valid_count)
    outraw = total/(vd+0.0001)
    if ll[5] != '-':
        out = list(outraw)
    else:
        out = list(outraw)[::-1]
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
    
def getsignal(inputfile,outputfile,BGmatrix,pcut,ncut,conserve,chipbw,mnasebw,dnase_span,mnase_span,gen,fetch_length=100):

    
 #   p=BwIO(pcut)
 #   chrom_len = {}
 #   for i in p.chromosomeTree['nodes']:
 #       chrom_len[i['key']] = i['chromSize']
    genome = twobitreader.TwoBitFile(gen)
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    conservebw = BigWigFile(open(conserve, 'rb'))
    chip_bw = BigWigFile(open(chipbw,'rb'))
    mnase_bw = BigWigFile(open(mnasebw,'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    pspan = dnase_span - ml/2
    inf.seek(0)
    pBG,nBG = readBG(BGmatrix)
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()
        if ll[0] == 'chrY':
            continue
        chrom = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        strand = ll[5]
        seq = genome[chrom][(start-pspan-3):(end + pspan+3)]
        pout = make_cut(pcutbw,ll,pspan,fetch_length)
        nout = make_cut(ncutbw,ll,pspan,fetch_length)
        conserveout = make_cut(conservebw,ll,pspan,fetch_length)
        chipout = make_chip(chip_bw,ll,len(conserveout))
        mnaseout = make_mnase(mnase_bw,ll,mnase_span)
        #print couserveout
        if strand == "-":
            pout,nout = nout,pout
        if pout == 'NA':
            continue        

        if 'N' in seq.upper():
            continue
        #print 1
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - 6):
            p.append(pBG[pseq[k:k+6].upper()])
            n.append(nBG[nseq[k:k+6].upper()])
        if strand != '-':
            pbglist = p
            nbglist = n
        else:
            pbglist = n[::-1]
            nbglist = p[::-1]
        newll = ll  + mnaseout + pout + nout + pbglist + nbglist + conserveout + chipout
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
                         
    optparser.add_option("-m","--mnasebw",dest="mnasebw",type="str",
                         help="")
    optparser.add_option("-c","--conserve",dest="conserve",type="str",default = "/mnt/Storage/home/huse/Data/conserv/hg19_v_conserv.bw",
                         help="default is /mnt/Storage/home/huse/Data/conserv/hg19_v_conserv.bw")
    optparser.add_option("--chipbw",dest="chipbw",type="str",
                         help="")
                         
    optparser.add_option("--Dspan",dest="Dspan",type="int",default = 25,
                         help="range of DNase signal,1bp per bin, default = 25 , means total length = 50")
    optparser.add_option("--Mspan",dest="Mspan",type="int",default = 1000,
                         help="range of MNase signal,10bp per bin, default = 1000 , means total length = 2000")

    optparser.add_option("--genome",dest="genome",type="str",default = '/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="2bit format")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    bgmatrix = options.bgmatrix
    pcut = options.pcut
    ncut = options.ncut
    genome = options.genome
    Dspan = options.Dspan
    Mspan = options.Mspan
    conserve = options.conserve
    chipbw = options.chipbw
    mnasebw = options.mnasebw

    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,bgmatrix,pcut,ncut,conserve,chipbw,mnasebw,Dspan,Mspan,genome,fetch_length=100)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
