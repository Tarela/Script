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

def sitepro_scan(pattern,peak,out,w_plus,w_minus,trunk,text):
    inf = open(pattern)
    pattern_plus = map(float,inf.readline().strip().split(","))
    pattern_minus = map(float,inf.readline().strip().split(","))
    all_sum = sum(pattern_plus)+sum(pattern_minus)
    p_plus = []
    p_minus= []
    for i in pattern_plus:
        p_plus.append(i/all_sum)
    for i in pattern_minus:
        p_minus.append(i/all_sum)
    inf.close()
    l = len(pattern_plus)
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
    footprint = []
#    count = 0
#    t=time.time()
    ls=[0]*44
    for line in inf:
 #       s=[]
        ll = line.split()
        if chrom_len1.has_key(ll[0])  and chrom_len2.has_key(ll[0]):
            p_sum = list(w_plus_H.summarize(ll[0],int(ll[1])-3-22,int(ll[2])+3+22,(int(ll[2])-int(ll[1])+6+44)).sum_data)
            m_sum = list(w_minus_H.summarize(ll[0],int(ll[1])-3-22,int(ll[2])+3+22,(int(ll[2])-int(ll[1])+6+44)).sum_data)
            last_start = "NA"
            last_end = "NA"
            last_value = "NA"
            for i in range(len(p_sum)-l):
                o_plus = map(int,p_sum[i:i+l])
                o_minus = map(int,m_sum[i:i+l])
                for k in range(len(o_plus)):
                    if o_plus[k] > trunk:
                        o_plus[k]=trunk
                    if o_minus[k] > trunk:
                        o_minus[k] = trunk
               # o_sum = sum(o_plus)+sum(o_minus)

                #print pattern_plus,p0
                score =  match_pattern(p_plus,p_minus,p0,p0,o_plus,o_minus,l)
    #            s.append(score)
                if i == 22:
                    footprint.append(ll+[score])
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
            #footprint.append(ll+[])
       # if count%100 ==0:
       #     print time.time()-t
       #     print ls
       #     t = time.time()
       # count += 1
    outf = open(out,'w')
    for fp in footprint:
        newline = "\t".join(map(str,fp))+"\n"
        outf.write(newline)
    outf.close()
    outf = open(text,'w')
    outf.write("\t".join(map(str,ls))+"\n")
    outf.close()
              #  print score



def Poisson_score(lam,k):
    return pow(lam,k)#/(reduce(lambda x,y:x*y, range(1,k+1)))

def match_pattern(pattern_plus,pattern_minus,pattern_plus0,pattern_minus0,observe_plus,observe_minus,l):
    ### l is the length of pattern
    observe_sum = sum(observe_plus)+sum(observe_minus)
    #print observe_sum
    lambda_plus = [x*observe_sum for x in pattern_plus]
    #print pattern_plus
    #print lambda_plus
    lambda_minus = [x*observe_sum for x in pattern_minus]
    lambda_plus0 = [x*observe_sum for x in pattern_plus0]
    lambda_minus0 = [x*observe_sum for x in pattern_minus0]
    #print lambda_plus
    #print lambda_minus
    #print lambda_plus0
    #print lambda_minus0
    f1=1
    f0=1
    for i in range(l):
        f1 = f1*Poisson_score(lambda_plus[i],observe_plus[i])*Poisson_score(lambda_minus[i],observe_minus[i])/(Poisson_score(lambda_plus0[i],observe_plus[i])*Poisson_score(lambda_minus0[i],observe_minus[i]))
    #print f1,f0
    #print math.log10(f1/f0)
      #  print f1,f0
    #print f1,f0
    if f1 == 0 :
        print f1
        return "NA"
    else:
        return math.log10(f1)


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
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/CFCE/Pair_end_DNase/Hiseq_pairend_DNase/Single_end_DNase_Data/mapping/25U_50U_single_plus_cutsite.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/CFCE/Pair_end_DNase/Hiseq_pairend_DNase/Single_end_DNase_Data/mapping/25U_50U_single_minus_cutsite.bw",
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

    if not pattern:
        optparser.print_help()
        sys.exit(1)
    
    sitepro_scan(pattern,interval,output,w_plus,w_minus,trunk,text)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

