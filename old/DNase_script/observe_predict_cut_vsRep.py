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
    
def getsignal(inputfile,outputfile,BGmatrix,pcutC,ncutC,pcut,ncut,pcut2,ncut2,chipbw,pspan,tspan,gen):

    
 #   p=BwIO(pcut)
 #   chrom_len = {}
 #   for i in p.chromosomeTree['nodes']:
 #       chrom_len[i['key']] = i['chromSize']
    genome = twobitreader.TwoBitFile(gen)
    pcutCbw = BigWigFile(open(pcutC, 'rb'))
    ncutCbw = BigWigFile(open(ncutC, 'rb'))
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    pcutbw2 = BigWigFile(open(pcut2, 'rb'))
    ncutbw2 = BigWigFile(open(ncut2, 'rb'))
    chip_bw = BigWigFile(open(chipbw,'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    pspan = pspan
    inf.seek(0)
    pBG,nBG = readBG(BGmatrix)
    outf = open(outputfile,'w')
    for line in inf:
        ll = line.split()

        chrom = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        strand = ll[5]#pspan = 5, tspan = 25
        
        seq = genome[chrom][(start-pspan-tspan-3):(end + pspan + tspan+3)]
        if start - pspan - tspan -3 < 0:
            continue
        if 'N' in seq.upper():
            continue
        pseq = seq[:-1]
        nseq = seq[1:] 

        poutC = list(pcutCbw.summarize(chrom,(start-pspan-tspan),(end+pspan+tspan),(end-start+2*pspan+2*tspan)).sum_data)
        pout = list(pcutbw.summarize(chrom,(start-pspan-tspan),(end+pspan+tspan),(end-start+2*pspan+2*tspan)).sum_data)
        pout2 = list(pcutbw2.summarize(chrom,(start-pspan-tspan),(end+pspan+tspan),(end-start+2*pspan+2*tspan)).sum_data)

        noutC = list(ncutCbw.summarize(chrom,(start-pspan-tspan),(end+pspan+tspan),(end-start+2*pspan+2*tspan)).sum_data)
        nout = list(ncutbw.summarize(chrom,(start-pspan-tspan),(end+pspan+tspan),(end-start+2*pspan+2*tspan)).sum_data)
        nout2 = list(ncutbw2.summarize(chrom,(start-pspan-tspan),(end+pspan+tspan),(end-start+2*pspan+2*tspan)).sum_data)

        pobsC = poutC[tspan:(tspan+end-start+2*pspan)]
        pobs = pout[tspan:(tspan+end-start+2*pspan)]
        pobs2 = pout2[tspan:(tspan+end-start+2*pspan)]

        nobsC = noutC[tspan:(tspan+end-start+2*pspan)]
        nobs = nout[tspan:(tspan+end-start+2*pspan)]
        nobs2 = nout2[tspan:(tspan+end-start+2*pspan)]
        pbias = []
        nbias = []
        for k in range(len(pseq)  +1 - 6):
            pbias.append(pBG[pseq[k:k+6].upper()])
            nbias.append(nBG[nseq[k:k+6].upper()])
        ppred = []
        npred = []
        for bp in range(len(pout)- 2*tspan):
        
            ptotal = sum(poutC[bp:(bp+2*tspan)])*1.0
            ntotal = sum(noutC[bp:(bp+2*tspan)])*1.0
        
            ppred.append(ptotal * pbias[bp+tspan]/sum(pbias[bp:(bp+2*tspan)]))
            npred.append(ntotal * nbias[bp+tspan]/sum(nbias[bp:(bp+2*tspan)]))

        if strand == "-":

            pobsC,nobsC = nobsC[::-1],pobsC[::-1]
            pobs,nobs = nobs[::-1],pobs[::-1]
            pobs2,nobs2 = nobs2[::-1],pobs2[::-1]
            ppred,npred = npred[::-1],ppred[::-1]
        
        chipout = make_chip(chip_bw,ll,1)[0]
        #print couserveout
      #  print tspan,(tspan+end-start+2*pspan),end,start,pspan
      #  print len(pobs),len(nobs),len(pobs2),len(nobs2),len(ppred),len(npred)
        newll = ll  + [chipout] + pobsC + nobsC + pobs + nobs + pobs2 + nobs2 + ppred + npred  
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
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_6mer_bg.txt",
                         help="")
    optparser.add_option("--positive_cut_combine",dest="pcutC",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_p.bw",
                         help="")
    optparser.add_option("--negative_cut_combine",dest="ncutC",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_m.bw",
                         help="")
    optparser.add_option("-p","--positive_cut",dest="pcut",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_r1_p.bw",
                         help="")
    optparser.add_option("-n","--negative_cut",dest="ncut",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_r1_m.bw",
                         help="")
    optparser.add_option("-q","--positive_2_cut",dest="p2cut",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_r2_p.bw",
                         help="")
    optparser.add_option("-l","--negative_2_cut",dest="n2cut",type="str",default="/mnt/Storage/home/huse/Project/DNase/Data/GM12878_UWDHS/GM12878_UWDHS_r2_m.bw",
                         help="")
    optparser.add_option("-c","--chipbw",dest="chipbw",type="str",
                         help="")
    optparser.add_option("-s","--pspan",dest="pspan",type="int",default = 5,
                         help="default = 5")
    optparser.add_option("-t","--tagcount_span",dest="tspan",type="int",default = 25,
                         help="default = 25 , means total length = 50")

    optparser.add_option("--genome",dest="genome",type="str",default = '/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="2bit format")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    bgmatrix = options.bgmatrix
    pcutC = options.pcutC
    ncutC = options.ncutC
    pcut = options.pcut
    ncut = options.ncut
    pcut2 = options.p2cut
    ncut2 = options.n2cut
    genome = options.genome
    tspan = options.tspan
    chipbw = options.chipbw
    pspan = options.pspan
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    getsignal(inputfile,outputfile,bgmatrix,pcutC,ncutC,pcut,ncut,pcut2,ncut2,chipbw,pspan,tspan,genome)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
