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
import math,time
from CistromeAP.taolib.CoreLib.BasicStat.Func import *
from CistromeAP.jianlib.BwReader import BwIO
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

def sitepro_scan(pattern,peak,out,w_plus,w_minus,trunk,text,w_chip):
    inf = open(pattern)
    pattern_plus_pmotif = map(float,inf.readline().strip().split(","))
    pattern_minus_pmotif = map(float,inf.readline().strip().split(","))
    pattern_plus_mmotif = map(float,inf.readline().strip().split(","))
    pattern_minus_mmotif = map(float,inf.readline().strip().split(","))

    all_sum_p = sum(pattern_plus_pmotif)+sum(pattern_minus_pmotif)
    all_sum_m = sum(pattern_plus_mmotif)+sum(pattern_minus_mmotif)

    p_plus = []
    p_minus= []
    m_plus = []
    m_minus = []
    for i in range(len(pattern_plus_pmotif)):
        p_plus.append(pattern_plus_pmotif[i]/all_sum_p)
        p_minus.append(pattern_minus_pmotif[i]/all_sum_p)
        m_plus.append(pattern_plus_mmotif[i]/all_sum_m)
        m_minus.append(pattern_minus_mmotif[i]/all_sum_m)
    inf.close()
    l = len(pattern_plus_pmotif)
    p0 = [1.0/(2*l)]*l
    inf = open(peak)
    p=BwIO(w_plus)
    q=BwIO(w_minus)
    chrom_len1 = {}
    chrom_len2 = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len1[i['key']] = i['chromSize']
    for i in q.chromosomeTree['nodes']:
        chrom_len2[i['key']] = i['chromSize']
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    w_chip_H=BigWigFile(open(w_chip, 'rb'))
    footprint = []
    ls=[0]*2*len(pattern_plus_pmotif)
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        if chrom_len1.has_key(ll[0])  and chrom_len2.has_key(ll[0]):
            DNase = float(w_plus_H.summarize(ll[0],(int(ll[1])+int(ll[2]))/2-50,(int(ll[1])+int(ll[2]))/2+50,1).sum_data) + float(w_minus_H.summarize(ll[0],(int(ll[1])+int(ll[2]))/2-50,(int(ll[1])+int(ll[2]))/2+50,1).sum_data) 
            Chip = float(w_chip_H.summarize(ll[0],int(ll[1]),int(ll[2]),1).sum_data) 
            p_sum = list(w_plus_H.summarize(ll[0],int(ll[1])-40-len(pattern_plus_pmotif),int(ll[1])+41+len(pattern_plus_pmotif),81+2*(len(pattern_plus_pmotif))).sum_data)
            m_sum = list(w_minus_H.summarize(ll[0],int(ll[1])-40-len(pattern_plus_pmotif),int(ll[2])+41+len(pattern_plus_pmotif),81+2*(len(pattern_plus_pmotif))).sum_data)
            last_start = "NA"
            last_end = "NA"
            last_value = "NA"
            for i in range(len(p_sum)-l):
                o_plus = map(int,p_sum[i:i+l])
                o_minus = map(int,m_sum[i:i+l])
#                for k in range(len(o_plus)):
#                    if o_plus[k] > trunk:
#                        o_plus[k]=trunk
#                    if o_minus[k] > trunk:
#                        o_minus[k] = trunk
                if ll[5]=="+":
                    score =  match_pattern(p_plus,p_minus,o_plus,o_minus,l,1)
                elif ll[5]=="-":
                    score =  match_pattern(m_plus,m_minus,o_plus,o_minus,l,1)
                if i == len(pattern_plus_pmotif):
                    footprint.append(ll+[score,DNase,Chip])
                if last_start == "NA" :
                    last_start = i
                    last_end = i+l
                    last_value = score
                elif score > last_value:
                    last_start = i
                    last_end = i+l
                    last_value = score
            
            if last_start ==0 and last_value ==0 :
                pass
            else:   
                ls[last_start]+=1
        
    outf = open(out,'w')
    for fp in footprint:
        newline = "\t".join(map(str,fp))+"\n"
        outf.write(newline)
    outf.close()
    outf = open(text,'w')
    outf.write("\t".join(map(str,ls))+"\n")
    outf.close()
              #  print score

def match_pattern(pattern_plus,pattern_minus,observe_plus,observe_minus,l,persudo):
    ### l is the length of pattern
    observe_sum = sum(observe_plus)+sum(observe_minus)
    if observe_sum == 0 :
        return "non"
    p0_sum = observe_sum + persudo*l
    lambda_plus = pattern_plus
    lambda_minus =  pattern_minus
    #print observe_plus
    lambda_plus0 = [(x*1.0+persudo)/p0_sum for x in observe_plus]
    lambda_minus0 = [(x*1.0+persudo)/p0_sum for x in observe_minus]
    f1=0
    #print lambda_plus
    for i in range(l):
    #    print lambda_plus[i]
    #    print lambda_minus[i]
        #print lambda_plus0[i]
        #print lambda_minus0[i]
      #  f1 = f1*Poisson_score(lambda_plus[i],observe_plus[i])*Poisson_score(lambda_minus[i],observe_minus[i])/(Poisson_score(lambda_plus0[i],observe_plus[i])*Poisson_score(lambda_minus0[i],observe_minus[i]))
        f1 = f1 + observe_plus[i]*(math.log10(lambda_plus[i])-math.log10(lambda_plus0[i])) + observe_minus[i]*(math.log10(lambda_minus[i])-math.log10(lambda_minus0[i]))
    return (f1/observe_sum)

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-p","--pattern",dest="pattern",type="str",
                         help="")
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("-t","--text",dest="text",type="str",
                         help="")
    optparser.add_option("-c","--chipseqbw",dest="chipseqbw",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_plus_50_100.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_minus_50_100.bw",
                         help="")
#========minor options=============
    optparser.add_option("--trunk",dest="trunk",type="float",default = 4.0,
                         help="")


    (options,args) = optparser.parse_args()

    pattern = options.pattern
    interval = options.interval
    output = options.output
    w_plus = options.w_plus
    w_minus = options.w_minus
    trunk = options.trunk
    text = options.text
    chipseqbw = options.chipseqbw
    if not pattern:
        optparser.print_help()
        sys.exit(1)
    
    sitepro_scan(pattern,interval,output,w_plus,w_minus,trunk,text,chipseqbw)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

