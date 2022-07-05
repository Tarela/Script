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

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def read_interval(peak):
    inf = open(peak)
    peak_dict = {}
    for line in inf:
        if not line.startswith("chr"):
            continue
        ll = line.split()
        if len(ll)<3:
            continue
    #    strand = 1
        #peakcenter = (int(ll[1])+int(ll[2]))/2
        if not peak_dict.has_key(ll[0]):
            peak_dict[ll[0]]=[]
        peak_dict[ll[0]].append((int(ll[1])+int(ll[2]))/2)
    for chrm in peak_dict.keys():
        peak_dict[chrm].sort()
    return peak_dict

def read_tss(Annotation):
    inf = open(Annotation)
    exp_dict = {}
    for line in inf:
        if not line.startswith("chr"):
            continue
        ll = line.split()
        if ll[5]=="+":
            tss = int(ll[1])
        elif ll[5]=="-":
            tss = int(ll[2])
        else:
            print "error"
            exit()
        if not exp_dict.has_key(ll[0]):
            exp_dict[ll[0]]=[]
        exp_dict[ll[0]].append(ll+[tss])
    for chrm in exp_dict.keys():
        exp_dict[chrm].sort(key = lambda x:x[6])
    return exp_dict

def distance(DHSdict,Expdict,out):
    outf = open(out,'w')
    
    for chrm in DHSdict.keys():
        Ds = DHSdict[chrm]
        if Expdict.has_key(chrm):
            Es = Expdict[chrm]
            Eindex=0
            Dindex=0
            mindis="NA"
            while 1:
                Rdis = ( Ds[Dindex][6]-Es[Eindex])
                dis = abs(Rdis)

                if Rdis >= 0:
                   # dis = abs( Ds[Dindex][6]-Es[Eindex])
                    if mindis == "NA" or mindis > dis:
                        mindis = dis
                    if Eindex < len(Es)-1:
                        Eindex += 1
                    else:
                        outf.write( "\t".join(Ds[Dindex][:6] + [str(mindis)]) + "\n")
                        if Dindex < len(Ds)-1:
                            Dindex += 1
                            mindis = "NA"
                            if Eindex > 0:
                                Eindex -= 1
                        else:
                            break
                else:
                    #dis = abs( Ds[Dindex][6] -Es[Eindex])
                    if mindis == "NA" or mindis > dis:
                        mindis = dis
                    outf.write( "\t".join(Ds[Dindex][:6] + [str(mindis)]) + "\n")
                    if Dindex < len(Ds)-1:
                        Dindex += 1
                        mindis="NA"
                        if Eindex > 0:
                            Eindex -= 1
                    else:
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
    optparser.add_option("-p","--peak",dest="peak",type="str",
                         help="")
    optparser.add_option("-g","--gene",dest="gene",type="str",
                         help="gene annotation in bed format")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="outputfile")
#========minor options=============

    (options,args) = optparser.parse_args()

    peak = options.peak
    gene = options.gene
    out = options.out

    if not peak or not gene:
        optparser.print_help()
        sys.exit(1)
    
    pd = read_interval(peak)
    gd = read_tss(gene)
    distance(gd,pd,out)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

