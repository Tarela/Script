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
    

def profile(inputfile,outputfile,mode,strand):
    A1=[0]*4##[ACGT]
    C1=[0]*4
    G1=[0]*4
    T1=[0]*4
    x = cis.interval(genome='hg19')
    x.chrom=[]
    x.start=[]
    x.end=[]
    inf = open(inputfile)
    for line in inf:
        ll = line.split()
        x.chrom.append(ll[0])
        if mode == "peak":
            x.start.append(int(ll[1]))
            x.end.append(int(ll[2])+2)
        elif mode =="cut":
            if strand == "+":
                x.start.append(int(ll[1]))
                x.end.append(int(ll[1])+2)
            else:
                x.start.append(int(ll[2]))
                x.end.append(int(ll[2])+2)
        else:
            print "mode wrong, only peak,cut availabe"
            exit()
    x.getSequence()
    for i in range(len(x.start)):
        s = string.upper(x.seq[i])
        if strand == "+":
            for j in range(len(s)-1):
                if s[j+1] == "A":
                    p1=0
                elif s[j+1] == "C":
                    p1=1
                elif s[j+1] == "G":
                    p1=2
                elif s[j+1] == "T":
                    p1=3
                else:
                    continue
                if s[j] == "A":
                    A1[p1]+=1
                elif s[j] == "C":
                    C1[p1]+=1
                elif s[j] == "G":
                    G1[p1]+=1
                elif s[j] == "T":
                    T1[p1]+=1
                else:
                    continue
        elif strand == "-":
            s = s[::-1]
            for j in range(len(s)-1):
                if s[j+1] == "A":
                    p1=3
                elif s[j+1] == "C":
                    p1=2
                elif s[j+1] == "G":
                    p1=1
                elif s[j+1] == "T":
                    p1=0
                else:
                    continue
                if s[j] == "A":
                    T1[p1]+=1
                elif s[j] == "C":
                    G1[p1]+=1
                elif s[j] == "G":
                    C1[p1]+=1
                elif s[j] == "T":
                    A1[p1]+=1
                else:
                    continue
        else:
            print "strand only + and - "
            exit()
    inf.close()
    outf = open(outputfile,'w')
#    outf.write("\t".join(['P','A','C','G','T'])+"\n")
#    outf.write("\t".join(map(str,['A']+A1))+"\n")
#    outf.write("\t".join(map(str,['C']+C1))+"\n")
#    outf.write("\t".join(map(str,['G']+G1))+"\n")
#    outf.write("\t".join(map(str,['T']+T1))+"\n")
    outf.write("\t".join(map(str,A1))+"\n")
    outf.write("\t".join(map(str,C1))+"\n")
    outf.write("\t".join(map(str,G1))+"\n")
    outf.write("\t".join(map(str,T1))+"\n")

    
#    outf.write("total\t"+str(sum(A1)+sum(C1)+sum(G1)+sum(T1))+"\n")
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
    optparser.add_option("-m","--mode",dest="mode",type="str",default="peak",
                         help="choose from peak and cut")
    optparser.add_option("-s","--strand",dest="strand",type="str",default="+",
                         help="choose from + and -")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outputfile = options.outputfile
    mode = options.mode
    strand = options.strand
    if not inputfile or not outputfile:
        optparser.print_help()
        sys.exit(1)
    
    profile(inputfile,outputfile,mode,strand)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

