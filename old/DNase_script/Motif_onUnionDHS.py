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

import os,sys,re,time
from optparse import OptionParser
import logging
import string

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def readDHS(DHSpeak):
    DHSdict = {}
    inf = open(DHSpeak)
    for line in inf:
        ll = line.split()
        if not DHSdict.has_key(ll[0]):
            DHSdict[ll[0]] = []
        DHSdict[ll[0]].append(ll)
    inf.close()
    for chrm in DHSdict.keys():
        DHSdict[chrm].sort(key = lambda x:int(x[1]))
    return DHSdict

def readmotif(motif):
    motifdict = {}
    inf = open(motif)
    for line in inf:
        ll = line.split()
        if not motifdict.has_key(ll[0]):
            motifdict[ll[0]] = []
        motifdict[ll[0]].append([int(ll[1]),int(ll[2])])
    inf.close()
    for chrm in motifdict.keys():
        motifdict[chrm].sort(key = lambda x:x[0])
    return motifdict

def overlap(DHSfile,motiffile,out):
    DHS = readDHS(DHSfile)
    motif = readmotif(motiffile)
    outf = open(out,'w')
    outf.write('chrm\tstart\tend\tunknown1\tunknown2\tstrand\t%s\n'%(motiffile.split('/')[-1].split('.')[0]).rstrip('_upper'))
    for chrm in sorted(DHS.keys()):
        if not motif.has_key(chrm):
            for dhs in DHS[chrm]:
                outf.write("\t".join(dhs+['0'])+"\n")
        else:
            Mindex=0
            Dindex=0
            dhs = DHS[chrm][Dindex] + [0]
            while 1:

                if int(DHS[chrm][Dindex][1]) > motif[chrm][Mindex][1] :
                    if Mindex < len(motif[chrm]) -1 :
                        Mindex += 1
                    else:
                        if Dindex < len(DHS[chrm]) - 1:
                            outf.write("\t".join(map(str,dhs))+"\n")
                            Dindex += 1
                            dhs = DHS[chrm][Dindex] + [0]
                        else:
                            outf.write("\t".join(map(str,dhs))+"\n")
                            break
                elif int(DHS[chrm][Dindex][1]) <= motif[chrm][Mindex][1] and int(DHS[chrm][Dindex][2]) >= motif[chrm][Mindex][1]:
                    dhs[-1] += 1
                    if Mindex < len(motif[chrm]) -1 :
                        Mindex += 1
                    else:
                        if Dindex < len(DHS[chrm]) - 1:
                            outf.write("\t".join(map(str,dhs))+"\n")
                            Dindex += 1
                            dhs = DHS[chrm][Dindex] + [0]
                        else:
                            outf.write("\t".join(map(str,dhs))+"\n")
                            break
                else:
                    if Dindex < len(DHS[chrm]) - 1:
                        outf.write("\t".join(map(str,dhs))+"\n")
                        Dindex += 1
                        dhs = DHS[chrm][Dindex] + [0]
                    else:
                        outf.write("\t".join(map(str,dhs))+"\n")
                        break
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
    optparser.add_option("-m","--motif",dest="motif",type="str",
                         help="")
    optparser.add_option("-p","--DHSpeak",dest="DHSpeak",type="str",default="/home/sh430/Data/DNase/ENCODEold/UnionDHS/hg19_intronintergenic_DHS.bed",
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    DHSpeak = options.DHSpeak
    motif = options.motif
    out = options.out

    if not motif:
        optparser.print_help()
        sys.exit(1)
    
    overlap(DHSpeak,motif,out)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

