'''
Created on XXXX-XX-XX

@author: Tarela
'''
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
import cistrome.regions as cis
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def readpattern(patternfile):
    p=[[],[],[],[]]
    inf = open(patternfile)
    a=inf.readline()
    p[0]=map(float,a.split())
    a=inf.readline()
    p[1]=map(float,a.split())
    a=inf.readline()
    p[2]=map(float,a.split())
    a=inf.readline()
    p[3]=map(float,a.split())    
    inf.close()
    return p

def profile(inputfile,outputfile,pattern,strand):
    x = cis.interval(genome='hg19')
    x.chrom=[]
    x.start=[]
    x.end=[]
    inf = open(inputfile)
    for line in inf:
        ll = line.split()
        x.chrom.append(ll[0])
        x.start.append(int(ll[1]))
        x.end.append(int(ll[2])+2)
    x.getSequence()
    scores = []
    for i in range(len(x.start)):
        score = []
        s = string.upper(x.seq[i])
        if strand == "+":
            for j in range(len(s)-2):
                if s[j] == "A":
                    p1=0
                elif s[j] == "C":
                    p1=1
                elif s[j] == "G":
                    p1=2
                elif s[j] == "T":
                    p1=3
                else:
                    break
                if s[j+1] == "A":
                    p2=0
                elif s[j+1] == "C":
                    p2=1
                elif s[j+1] == "G":
                    p2=2
                elif s[j+1] == "T":
                    p2=3
                else:
                    break
                prob = pattern[p1][p2]
                score.append(prob)
        elif strand == "-":
            for j in range(1,len(s)-1):
                if s[j+1]=="A":
                    p1=3
                elif s[j+1] == "C":
                    p1=2
                elif s[j+1] == "G":
                    p1=1
                elif s[j+1] == "T":
                    p1=0
                else:
                    break
                if s[j] == "A":
                    p2=3
                elif s[j] == "C":
                    p2=2
                elif s[j] == "G":
                    p2=1
                elif s[j] == "T":
                    p2=0
                else:
                    break
                prob = pattern[p1][p2]
                score.append(prob)
        scores.append(score)
    inf.close()
    outf = open(outputfile,'w')
    for score in scores:
        outf.write("\t".join(map(str,score))+"\n")
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
    optparser.add_option("-s","--strand",dest="strand",type="str",default="+",
                         help="choose from + and -")
    optparser.add_option("-p","--pattern",dest="patternfile",type="str",
                         help="file for the matrix")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    strand = options.strand
    patternfile = options.patternfile
    if not inputfile or not outputfile:
        optparser.print_help()
        sys.exit(1)
    pattern=readpattern(patternfile)
    profile(inputfile,outputfile,pattern,strand)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

