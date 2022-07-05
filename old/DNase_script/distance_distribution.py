'''
Created on 

@author: 
'''
#!/usr/bin/env python
#Time-stamp:<>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import string
import time
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def read_interval(peak,span):
    inf = open(peak)
    peak_dict = {}
    for line in inf:
        if not line.startswith("chr"):
            continue
        ll = line.split()
        if len(ll)<3:
            continue
        if not peak_dict.has_key(ll[0]):
            peak_dict[ll[0]]=[[],[],[],[]]
        peak_dict[ll[0]][0].append((int(ll[1])+int(ll[2]))/2-span)
        peak_dict[ll[0]][1].append((int(ll[1])+int(ll[2]))/2+span)
        peak_dict[ll[0]][2].append((int(ll[1])+int(ll[2]))/2)
        if len(ll) >= 6 and ll[5]=="-":
            peak_dict[ll[0]][3].append("-")
        else:
            peak_dict[ll[0]][3].append("+")
    return peak_dict

def read_tag(tag):
    inf = open(tag)
    tag_dict={}
    for line in inf:
        if not line.startswith("chr"):
            continue
        ll = line.split()
        if len(ll)<3:
            print "your bed less than 3 column ~ seems wrong file"
            exit()
        if not tag_dict.has_key(ll[0]): ## do not have chrom before
            tag_dict[ll[0]]=[]
        tag_dict[ll[0]].append((int(ll[2])+int(ll[1]))/2)
    return tag_dict


def getoverlap(peaklist,tagiter,current_result):
    i=0
    pall=len(peaklist)
    t=tagiter.next()
 ##peak : [[left],[right],[mid],[strand]]
 ##tag  : [mid] 
    while 1:
        if t < peaklist[i][0]:
            try:
                t=tagiter.next()
                resT=0
            except StopIteration:
                break
        elif t >= peaklist[i][0] and t <= peaklist[i][1]:
            if peaklist[i][3]=="+":
                current_result[t-peaklist[i][2]] += 1
            else:
                current_result[peaklist[i][2]-t] += 1
            ii = i+1
            if ii < pall:
                while t >= peaklist[ii][0] and t <= peaklist[ii][1]:
                    if peaklist[ii][3]=="+":
                        current_result[t-peaklist[ii][2]] += 1
                    else:
                        current_result[peaklist[ii][2]-t] += 1
                    if ii <  pall -1 :
                        ii += 1
                    else:   
                        break
            try:
                t=tagiter.next()
            except StopIteration:
                if i < pall - 1:
                    i += 1
                else:
                    break
        else:
            if i < pall-1:
                i+=1
                resP=0
            else:
                break
    return current_result

def cal_dis(peaklist,tgl,current_result):
##peak = motif : [[mid],[strand]]
##tag = footprint  : [mid]   
    peak_iter = (sorted(zip(peaklist[0],peaklist[1],peaklist[2],peaklist[3]),key=lambda x:x[0]))
    tag_iter = iter(sorted(tgl))
    current_result=getoverlap(peak_iter,tag_iter,current_result)
    return current_result

def distance_length(peak_dict,tag_dict,span):
    result_dict = {}
    for i in range(-span,span+1):
        result_dict[i]=0
    for chrom in peak_dict.keys():  ##peak : [[left],[right],[mid],[strand]]
                                    ##tag  : [mid]                
        if not tag_dict.has_key(chrom):
            continue
        result_dict = cal_dis(peak_dict[chrom],tag_dict[chrom],result_dict)
    return result_dict

def output(result_dict,out):
    outf = open(out,'w')
    x=[]
    y=[]
    for dis in sorted(result_dict.keys()):
        x.append(dis)
        y.append(result_dict[dis])
    lx = map(str,x)
    ly = map(str,y)
    linex = "\t".join(lx)+"\n"
    liney = "\t".join(ly)+"\n"
    outf.write(linex)
    outf.write(liney)
    outf.close()
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-p interval> <-t pair_end tag> <-o output_name> [options]"
    description = """V plot"""

    optparser = OptionParser(version="%prog >.<",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-p","--motif_interval",dest="peak",type="str",
                         help="your interval for xaxis of the plot , the script takes the center of each region for caculating the distance")
    optparser.add_option("-t","--Footprint",dest="tag",type="str",
                         help="tag file in 6 column bed format ")
    optparser.add_option("-o","--output_name",dest="out",type="str",
                         help="name of your outputfile,2 line , x and y for R plot")
    optparser.add_option("--span",dest="span",type="int",default=100,
                         help="length of each region from the region center")

    (options,args) = optparser.parse_args()

    peak = options.peak
    tag = options.tag
    out = options.out
    span = abs(options.span)
    if not tag or not peak:
        optparser.print_help()
        sys.exit(1)
    c=time.time()
    tag_dict = read_tag(tag)
    print "read tag ",time.time()-c
    c=time.time()
    b=time.time()
    peak_dict = read_interval(peak,span)
    print "read interval ",time.time()-b
    e=time.time()
    result_dict = distance_length(peak_dict,tag_dict,span)
    print "caculate ",time.time()-e
    # print result_dict[0]
    output(result_dict,out)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

