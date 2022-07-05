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
    
def getsignal(inputfile,outputfile,atacBGmatrix,dnaseBGmatrix,atacpcut,atacncut,atacshortpcut,atacshortncut,dnasepcut,dnasencut,conserve,dnase_span,gen,fetch_length=100):

    
 #   p=BwIO(pcut)
 #   chrom_len = {}
 #   for i in p.chromosomeTree['nodes']:
 #       chrom_len[i['key']] = i['chromSize']
    genome = twobitreader.TwoBitFile(gen)
    atacpcutbw = BigWigFile(open(atacpcut, 'rb'))
    atacncutbw = BigWigFile(open(atacncut, 'rb'))
    atacshortpcutbw = BigWigFile(open(atacshortpcut, 'rb'))
    atacshortncutbw = BigWigFile(open(atacshortncut, 'rb'))
    dnasepcutbw = BigWigFile(open(dnasepcut, 'rb'))
    dnasencutbw = BigWigFile(open(dnasencut, 'rb'))

    conservebw = BigWigFile(open(conserve, 'rb'))
    inf = open(inputfile)    
    testll = inf.readline().split()
    ml = int(testll[2]) - int(testll[1])
    pspan = dnase_span - ml/2
    inf.seek(0)
    atacpBG,atacnBG = readBG(atacBGmatrix)
    dnasepBG,dnasenBG = readBG(dnaseBGmatrix)
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
        atacpout = make_cut(atacpcutbw,ll,pspan,fetch_length)
        atacnout = make_cut(atacncutbw,ll,pspan,fetch_length)
        atacshortpout = make_cut(atacshortpcutbw,ll,pspan,fetch_length)
        atacshortnout = make_cut(atacshortncutbw,ll,pspan,fetch_length)
        dnasepout = make_cut(dnasepcutbw,ll,pspan,fetch_length)
        dnasenout = make_cut(dnasencutbw,ll,pspan,fetch_length)
        conserveout = make_cut(conservebw,ll,pspan,fetch_length)
        #print couserveout
        if strand == "-":
            atacpout,atacnout = atacnout,atacpout
            atacshortpout,atacshortnout = atacshortnout,atacshortpout
            dnasepout,dnasenout = dnasenout,dnasepout
            
        if atacpout == 'NA':
            continue        

        if 'N' in seq.upper():
            continue
        #print 1
        pseq = seq[:-1]
        nseq = seq[1:]
        atacp=[]
        atacn=[]
        dnasep = []
        dnasen = []
        for k in range(len(pseq)  +1 - 6):
            atacp.append(atacpBG[pseq[k:k+6].upper()])
            atacn.append(atacnBG[nseq[k:k+6].upper()])
            dnasep.append(dnasepBG[pseq[k:k+6].upper()])
            dnasen.append(dnasenBG[nseq[k:k+6].upper()])
        if strand != '-':
            atacpbglist = atacp
            atacnbglist = atacn
            dnasepbglist = dnasep
            dnasenbglist = dnasen
        else:
            atacpbglist = atacn[::-1]
            atacnbglist = atacp[::-1]
            dnasepbglist = dnasen[::-1]
            dnasenbglist = dnasep[::-1]
        newll = ll  + atacpout + atacnout + atacshortpout + atacshortnout + atacpbglist + atacnbglist + dnasepout + dnasenout + dnasepbglist + dnasenbglist + conserveout 
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
    optparser.add_option("--atacbgmatrix",dest="atacBGmatrix",type="str",
                         help="")
    optparser.add_option("--dnasebgmatrix",dest="dnaseBGmatrix",type="str",
                         help="")
    optparser.add_option("--atacpcut",dest="atacpcut",type="str",
                         help="")
    optparser.add_option("--atacncut",dest="atacncut",type="str",
                         help="")
    optparser.add_option("--atacshortpcut",dest="atacshortpcut",type="str",
                         help="")
    optparser.add_option("--atacshortncut",dest="atacshortncut",type="str",
                         help="")
    optparser.add_option("--dnasepcut",dest="dnasepcut",type="str",
                         help="")
    optparser.add_option("--dnasencut",dest="dnasencut",type="str",
                         help="")

    optparser.add_option("-c","--conserve",dest="conserve",type="str",default = "/mnt/Storage/home/huse/Data/conserv/hg19_v_conserv.bw",
                         help="default is /mnt/Storage/home/huse/Data/conserv/hg19_v_conserv.bw")
    optparser.add_option("--Dspan",dest="Dspan",type="int",default = 50,
                         help="range of DNase signal,1bp per bin, default = 50 , means total length = 100")

    optparser.add_option("--genome",dest="genome",type="str",default = '/mnt/Storage/home/huse/Data/Genome/hg19/hg19.2bit',
                         help="2bit format")


#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    genome = options.genome
    Dspan = options.Dspan
    conserve = options.conserve

    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    getsignal(inputfile,outputfile,options.atacBGmatrix,options.dnaseBGmatrix,options.atacpcut,options.atacncut,options.atacshortpcut,options.atacshortncut,options.dnasepcut,options.dnasencut,conserve,Dspan,genome,fetch_length=100)    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
