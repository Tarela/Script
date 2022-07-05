#!/usr/bin/env python
'''
Created on XXXX-XX-XX

@author: Tarela
'''
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging,time
import string
import twobitreader

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def CpGrate(SEQ):
    Ccount = 0
    Gcount = 0
    CpGcount = 0
    for i in SEQ:
        if i.upper()=="G":# or i.upper()=="C":
            Gcount +=1
        elif i.upper() == "C":
            Ccount += 1
    for i in range(len(SEQ)-1):
        bp = SEQ[i:(i+2)]
        if bp.upper() == "CG":
            CpGcount += 1
    if Ccount*Gcount == 0:
        rate = 0
#        print SEQ
    else:
        rate = CpGcount *1.0 * len(SEQ) / (Ccount*Gcount)
    return rate
def CpG(fullSEQ):
    Mrate = []
    for i in range((len(fullSEQ)-500)/50 + 1):
        seq = fullSEQ[(i*50):(i*50+500)]
        if not len(seq) == 500: 
            print seq
        Mrate.append(CpGrate(seq))
    return max(Mrate)
def CGcontent(fullSEQ):
    GC=0
    for i in fullSEQ:
        if i.upper()=="G" or i.upper()=="C":
                GC +=1
    return str(GC*1.0/(len(fullSEQ)))


def get_CpG(inputfile,output,extend,sequence):
    genome = twobitreader.TwoBitFile(sequence) 
#print genome['chrM'][0:(dict['chrM']-2)]

    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        center = (int(ll[1])+int(ll[2]))/2
        seq1 = genome[ll[0]][max(0,center-extend):(center+extend)] #  (bwHandle.summarize(ll[0],max(0,center-extup),center+extdown,n))
        #signal=(bwHandle.summarize(ll[0],max(0,center-extend),center+extend,n))
        newll = ll+[CpG(seq1),CpGrate(seq1),CGcontent(seq1)]
        #print ll,len(seq)
        #print len(seq), ll[0],  center-extup, center+extdown
        #newll = ll +[R]
        #if ll[5] == "-":
        #    newll = ll+R[::-1]
        #else:
        #    newll = ll+R
        outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()



# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """add sequence feature following each line,(ext<-center->ext), newll  = ll+[localMaxCpG,aveCpG,GCcontent]"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("-g","--genome",dest="genome",type="str",
                         help="2bit genome sequence")
    optparser.add_option("--ext",dest="EXT",type="int",default = 1000,
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    sequence = options.genome
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    get_CpG(inputfile,output,options.EXT,sequence)
    #get_signal(inputfile,output,signalbw,extdown,extup,number)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
