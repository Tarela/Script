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
    
    
def getsignal(inputfile,outputfile,cut,pcut,ncut,pspan,tspan,BGmatrix,gen,left,right,pseudo):

    genome = twobitreader.TwoBitFile(gen)
#    p=BwIO(pcut)
#    chrom_len = {}
#    for i in p.chromosomeTree['nodes']:
#        chrom_len[i['key']] = i['chromSize']
    cutbw = BigWigFile(open(cut, 'rb'))
    pcutbw = BigWigFile(open(pcut, 'rb'))
    ncutbw = BigWigFile(open(ncut, 'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    inf.seek(0)
    outf = open(outputfile,'w')
    pBG,nBG = readBG(BGmatrix)
    for line in inf:
        ll = line.split()

        chrom = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        strand = ll[5]
        
        ## make tag count
        cut = list(cutbw.summarize(ll[0],int(ll[1]) + ml/2 -pspan ,int(ll[1]) + ml/2 +pspan ,2*pspan).sum_data)
        TC = sum(cut)
        ## make L, C , R
        pout = list(pcutbw.summarize(chrom, start - ml - tspan , end + ml + tspan , 3*ml + 2*tspan ).sum_data)
        nout = list(ncutbw.summarize(chrom, start - ml - tspan , end + ml + tspan , 3*ml + 2*tspan ).sum_data)
        
        if strand != '-':
            pass
        else:
            pout,nout = nout[::-1],pout[::-1]

      #  L = pcut[:ml] + ncut[:ml]
      #  C = pcut[ml: (ml*2)] + ncut[ml: (ml*2)]
      #  R = pcut[(ml*2):(ml*3)] + ncut[(ml*2):(ml*3)]
        ## make bias for L,C,R
        seq = genome[ll[0]][(int(ll[1]) - ml - tspan - left):(int(ll[2]) + ml + tspan + right)]
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
        p_predict = []
        n_predict = []
#        p_uniform = []
#        n_uniform = []
        p_LR = []
        n_LR = []
        
        for i in range(len(pout)-2*tspan):
            p_ob = pout[i+tspan]
            n_ob = nout[i+tspan]
            ptotal = sum(pout[i:(i+tspan*2)])
            ntotal = sum(nout[i:(i+tspan*2)])
            pbias = pbglist[i+tspan]
            nbias = nbglist[i+tspan]
            pbias_total = sum(pbglist[i:(i+tspan*2)])
            nbias_total = sum(nbglist[i:(i+tspan*2)])
            pbias_por = pbias*1.0/pbias_total
            nbias_por = nbias*1.0/nbias_total
            ppredict = pbias_por * ptotal
            npredict = nbias_por * ntotal
#            puniform = ptotal*1.0/(2*tspan)
#            nuniform = ntotal*1.0/(2*tspan)

            pLR = math.log(p_ob+pseudo) - math.log(ppredict+pseudo)
            nLR = math.log(n_ob+pseudo) - math.log(npredict+pseudo)

            p_cleavage.append(p_ob)
            n_cleavage.append(n_ob)
            p_bias.append(pbias)
            n_bias.append(nbias)
            p_predict.append(ppredict)
            n_predict.append(npredict)
#            p_uniform.append(puniform)
#            n_uniform.append(nuniform)
            p_LR.append(pLR)    
            n_LR.append(nLR)

        Lob = p_cleavage[:ml] + n_cleavage[:ml]
        Cob = p_cleavage[ml:(ml*2)] + n_cleavage[ml:(ml*2)]
        Rob = p_cleavage[(ml*2):(ml*3)] + n_cleavage[(ml*2):(ml*3)]
        
        Llr = p_LR[:ml] + n_LR[:ml]
        Clr = p_LR[ml:(ml*2)] + n_LR[ml:(ml*2)]
        Rlr = p_LR[(ml*2):(ml*3)] + n_LR[(ml*2):(ml*3)]

        FOS = -1*( (sum(Cob)+1)/(sum(Rob)+1) + (sum(Cob)+1)/(sum(Lob)+1) ) 
        lr_score = sum(Llr)+sum(Rlr)-sum(Clr)
        
        newll = ll + [TC,FOS,lr_score]
        outf.write("\t".join(map(str,newll))+"\n")

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
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outputfile",dest="outputfile",type="str",
                         help="")                         

    optparser.add_option("-g","--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="")

    optparser.add_option("-w","--cutbw",dest="cut",type="str",
                         help="")
    optparser.add_option("-p","--pcutbw",dest="pcut",type="str",
                         help="")
    optparser.add_option("-n","--ncutbw",dest="ncut",type="str",
                         help="")                         
    optparser.add_option("-s","--plot_span",dest="pspan",type="int",default = 100,
                         help="window size for tagcount and dymDHS, default = 100 , means total length = 200")
    optparser.add_option("-t","--count_span",dest="tspan",type="int",default = 25,
                         help="window size for window count to get predicted cut, default = 25 , means total length = 50")
    optparser.add_option("--left",dest="left",type="int",default = 3,
                         help="")
    optparser.add_option("--right",dest="right",type="int",default =3,
                         help="")                         
    optparser.add_option("--pseudo",dest="pseudo",type="int",default =1,
                         help="")                         

 
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    cutbw = options.cut
    pcutbw = options.pcut
    ncutbw = options.ncut
    BGmatrix = options.bgmatrix
    pspan = options.pspan
    tspan = options.tspan
    gen = options.genome
    left = options.left
    right = options.right
    pseudo = options.pseudo
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,cutbw,pcutbw,ncutbw,pspan,tspan,BGmatrix,gen,left,right,pseudo)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
