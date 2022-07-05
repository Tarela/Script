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
import copy,time,random,numpy

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def maek_sti(cleavage,seqbias):
    c_sum = sum(cleavage)
    s_sum = sum(seqbias)
    stimulate_seqbias = [0]*len(seqbias)
    proportion_seqbias = numpy.array(seqbias)/s_sum
    cumsum_seqbias = proportion_seqbias.cumsum()
    #print seqbias
   # print proportion_seqbias
   # print cumsum_seqbias
    for i in range(c_sum):
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


def stimulate(inputdata,out):#,cal_length):
    inf = open(inputdata)
    outf = open(out,'w')
    for line in inf:
        ll = line.split()
        ### read data
        motif = ll[:8]
        cal_length = (len(ll)-8)/6
        positive_cut = map(int,map(float,ll[8:(8+cal_length)]))
        negative_cut = map(int,map(float,ll[(8+cal_length):(8+cal_length*2)]))
        positive_bias = map(float,ll[(8+cal_length*4):(8+cal_length*5)])
        negative_bias = map(float,ll[(8+cal_length*5):(8+cal_length*6)])
        positive_sti = maek_sti(positive_cut,positive_bias)
        negative_sti = maek_sti(negative_cut,negative_bias)
        #print positive_sti
        #print positive_bias
        #print motif
        newll = ll + positive_sti + negative_sti
        outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
    inf.close()
        ### get sum , strand separately
        #positive_cut_sum = sum(positive_cut)
        #negative_cut_sum = sum(negative_cut)
        #positive_bias_sum = sum(positive_bias)
        #negative_bias_sum = sum(negative_bias)
        
    

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputdata",dest="inputdata",type="str",
                         help="MDC cleavage + seqbias data , generate by (daisy) Scan_matrix.py")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="output file , seqbias part -> stimulate")
#    optparser.add_option("-l","--cal_length",dest="cal_length",type="int",default=50,
#                         help="region length for calculate , AR is 50 , CTCF is 100")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    inputdata = options.inputdata
    out = options.out
  #  cal_length = options.cal_length
    if not inputdata:
        optparser.print_help()
        sys.exit(1)
    
    stimulate(inputdata,out)#,cal_length)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


