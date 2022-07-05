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

def read_interval(peak,span,orient):
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
            peak_dict[ll[0]]=[[],[],[],[]]
        peak_dict[ll[0]][0].append((int(ll[1])+int(ll[2]))/2-span)
        peak_dict[ll[0]][1].append((int(ll[1])+int(ll[2]))/2+span)
        peak_dict[ll[0]][2].append((int(ll[1])+int(ll[2]))/2)
        if orient and len(ll) >= 6 and ll[5]=="-":
            peak_dict[ll[0]][3].append("-")
        else:
            peak_dict[ll[0]][3].append("+")

    return peak_dict


def read_tag(tag,maxlen,minlen):
    inf = open(tag)
   # tag_dict = {}#{[pos start],[pos end],[neg start],[neg end]}
    tag_dict={}
    #neg_start={}
    #neg_end={}
    a=time.time()
    for line in inf:
        if not line.startswith("chr"):
            continue
        ll = line.split()
        if len(ll)<3:
            print "your bed less than 3 column ~ seems wrong file"
            exit()
        taglen = int(ll[2])-int(ll[1])
        if taglen > maxlen or taglen < minlen:
            continue
        if not tag_dict.has_key(ll[0]): ## do not have chrom before
            tag_dict[ll[0]]=[[],[],[],[],[]]
        tag_dict[ll[0]][0].append(int(ll[1]))
        tag_dict[ll[0]][1].append(int(ll[2]))
        tag_dict[ll[0]][2].append(int(ll[2])-int(ll[1]))
        tag_dict[ll[0]][3].append((int(ll[2])+int(ll[1]))/2)
        if len(ll)>=6 and ll[5] in ["+","-"]:
            tag_dict[ll[0]][4].append(ll[5])
        else:
            tag_dict[ll[0]][4].append("+")
    return tag_dict

def search(p, search_list):### number and search_list should be int and list of int,nearest dis is the caculate range
    maxi=(len(search_list)-1)
    mini=0
    while search_list[((maxi+mini)/2)][0] < p[0] or search_list[((maxi+mini)/2)][0] > p[1]:
        if search_list[((maxi+mini)/2)][0] > p[1]:
            maxi = (maxi+mini)/2
        else:
            mini = (maxi+mini)/2
        if maxi-mini <= 1:
            if search_list[maxi][0] < p[0] or search_list[maxi][0] > p[1]:
                if search_list[mini][0] < p[0] or search_list[mini][0] > p[1]:
                    return "NA"
                else:
                    return mini
            else:
                return maxi
    return (maxi+mini)/2

def getoverlap(peaklist,tagiter,current_result):
#    p_res = []
#    t_res = []
#    p=peakiter.next()
    i=0
    pall=len(peaklist)
    t=tagiter.next()
#    resP=0
#    resT=0
 #   overP=0
 ##peak : [[left],[right],[mid],[strand]]
 ##tag  : [[left],[right],[length],[mid],[strand]] 
 #   overT=0
    while 1:
 #       print i,peaklist[i]
        if t[0] < peaklist[i][0]:
            try:
                t=tagiter.next()
                resT=0
            except StopIteration:
                break
        elif t[0] >= peaklist[i][0] and t[0] <= peaklist[i][1]:
            if peaklist[i][3]=="+":
                current_result[t[0]-peaklist[i][2]][t[1]] += 1
            else:
                current_result[peaklist[i][2]-t[0]][t[1]] += 1
            ii = i+1
            if ii < pall:
                while t[0] >= peaklist[ii][0] and t[0] <= peaklist[ii][1]:
                    if peaklist[ii][3]=="+":
                        current_result[t[0]-peaklist[ii][2]][t[1]] += 1
                    else:
                        current_result[peaklist[ii][2]-t[0]][t[1]] += 1
                    if ii <  pall -1 :
                        ii += 1
                    else:   
                        break
#            if resT ==0:
#                t_res.append(t)
#                resT=1
#            if resP ==0:
#                p_res.append(peaklist[i])
#                resP=1
            try:
                t=tagiter.next()
               # resT=0
            except StopIteration:
                #try:
                    #p=peakiter.next()
                if i < pall - 1:
                    i += 1
               #     resP=0
                #except StopIteration:
                else:
                    break
        else:
            #try :
            #    p=peakiter.next()
            #    resP=0
            #except StopIteration:
            #    break
            if i < pall-1:
                i+=1
                resP=0
            else:
                break
    return current_result
def getoverlap_iter(peakiter,tagiter):
    p_res = []
    t_res = []
    p=peakiter.next()
#    i=0
#    pall=len(peaklist)
    t=tagiter.next()
    resP=0
    resT=0
 #   overP=0
 #   overT=0
    while 1:
        if t[0] < p[0]:
            try:
                t=tagiter.next()
                resT=0
            except StopIteration:
                break
        elif t[0] >= p[0] and t[0] <= p[1]:
            if resT ==0:
                t_res.append(t)
                resT=1
            if resP ==0:
                p_res.append(p)
                resP=1
            try:
                t=tagiter.next()
                resT=0
            except StopIteration:
                try:
                    p=peakiter.next()
                #if i < pall:
                #    i += 1
                    resP=0
                except StopIteration:
                #else:
                    break
        else:
            try :
                p=peakiter.next()
                resP=0
            except StopIteration:
                break
            #if i<pall:
            #    i+=1
            #    resP=0
            #else:
            #    break
    return p_res,t_res
def cal_dis(peaklist,tgl,current_result,tagpos):
##peak : [[left],[right],[mid],[strand]]
##tag  : [[left],[right],[length],[mid],[strand]]   
    peak_iter = (sorted(zip(peaklist[0],peaklist[1],peaklist[2],peaklist[3]),key=lambda x:x[0]))
    taglist = zip(tgl[0],tgl[1],tgl[2],tgl[3],tgl[4])
    if tagpos == "mid":
        tag_iter = iter(sorted([[w[3],w[2]] for w in taglist],key=lambda x:x[0]))
    elif tagpos == "5":
        tag_iter = iter(sorted([[w[0],w[2]] for w in taglist if w[4]=="+"]+[[w[1],w[2]] for w in taglist if w[4]=="-"],key=lambda x:x[0]))
    elif tagpos == "3":
        tag_iter = iter(sorted([[w[1],w[2]] for w in taglist if w[4]=="+"]+[[w[0],w[2]] for w in taglist if w[4]=="-"],key=lambda x:x[0]))
    gg=time.time()
    current_result=getoverlap(peak_iter,tag_iter,current_result)
    print "overlap time ",time.time()-gg
#    if len(p_res) ==0 or len(t_res)==0:
#        return current_result
#    p_res_iter = iter(p_res)
#    tag_number = len(t_res)
#    print tag_number
#    print len(p_res)
#    while 1:
#        try :
#            p = p_res_iter.next()
#        except StopIteration:
#            break
#        hit_index = search(p,t_res) 
    #    print p
    #    print t_res
#        if hit_index == "NA":
#            pass
#        else:
   #         print hit_index
#            for i in range(hit_index,tag_number):
#                if t_res[i][0] <= p[1] and t_res[i][0] >= p[0]:
#                    if p[3]=="-":
#                        current_result[p[2]-t_res[i][0]][t_res[i][1]]+=1
#                    else:
#                        current_result[t_res[i][0]-p[2]][t_res[i][1]]+=1
#                else:
#                    break
#            for j in range(hit_index,0,-1):
#                if  t_res[j][0] <= p[1] and t_res[j][0] >= p[0]:
#                    if p[3]=="-":
#                        current_result[p[2]-t_res[j][0]][t_res[j][1]]+=1
#                    else:
#                        current_result[t_res[j][0]-p[2]][t_res[j][1]]+=1
#                else:
#                    break 
    return current_result

def distance_length(peak_dict,tag_dict,tagpos,span,maxlen,minlen):
    result_dict = {}
    #taglenrange = maxlen-minlen
    #if judge_zero(ns):
    #    tagpos == "mid"
    for i in range(-span,span+1):
        result_dict[i]=[0]*(maxlen+1)
#    print result_dict[0]
    for chrom in peak_dict.keys():  ##peak : [[left],[right],[mid],[strand]]
                                    ##tag  : [[left],[right],[length],[mid],[strand]]                
        if not tag_dict.has_key(chrom):
            continue
        #plot_tag_list = []
        result_dict = cal_dis(peak_dict[chrom],tag_dict[chrom],result_dict,tagpos)
        #result_dict = cal_dis(plus_peak,nloci,nlength,result_dict,span,"+")    
        #result_dict = cal_dis(minus_peak,ploci,plength,result_dict,span,"-")
        #result_dict = cal_dis(minus_peak,nloci,nlength,result_dict,span,"-")
    return result_dict

def output(result_dict,out,minlen,maxlen,span):
    outf = open(out,'w')
    for taglen in range(minlen,maxlen+1):
        newll = []
        for pos in range(-span,span+1):
            newll.append(result_dict[pos][taglen])
      #      if result_dict[pos][taglen] != 0:
      #          print result_dict[pos][taglen]
        newll = map(str,newll)
        newline = "\t".join(newll)+"\n"
        outf.write(newline)
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
    optparser.add_option("-p","--peak_interval",dest="peak",type="str",
                         help="your interval for xaxis of the plot , the script takes the center of each region for caculating the distance")
    optparser.add_option("-t","--tag_pairend",dest="tag",type="str",
                         help="tag file in 6 column bed format ")
    optparser.add_option("-o","--output_name",dest="out",type="str",
                         help="name of your outputfile, the output file is a matrix, whose row for tag length and column for distance from center of each interval, the value in each cell represent the tag count in given taglength&distance")
    optparser.add_option("--span",dest="span",type="int",default=250,
                         help="length of each region from the region center , keep it same as --maxlen")
    optparser.add_option("--orient-peak",action="store_true",dest="orient",
                         help="If set, the direction (+/-) is considered in profiling. If no strand info given in the BED, this option is ignored. default:True")    
    optparser.add_option("--tag-pos",dest="tagpos",type="str",default = "mid",
                         help="""restriction enzyme end, the position of tag that caculate the distance from the center of interval , choice : [ 5 : 5\' of tag(only when tag has strand information) , 3 : 3\' of tag(only when tag has strand information) , mid : middle of tag , default : mid ]""")
    optparser.add_option("--maxlen",dest="maxlen",type="int",default=250,
                         help="the maximum of the distance between each tag & region")
    optparser.add_option("--minlen",dest="minlen",type="int",default=0,
                         help="the minimum of the distance between each tag & region")

    (options,args) = optparser.parse_args()

    peak = options.peak
    tag = options.tag
    out = options.out
    span = abs(options.span)
    orient = options.orient
    tagpos = options.tagpos
    maxlen = options.maxlen
    minlen = options.minlen
    if not tag or not peak:
        optparser.print_help()
        sys.exit(1)
    if not tagpos in ["5","3","mid","all"]:
        print """tagpos selection only in ["5","3","mid","all"]"""
        optparser.print_help()
        sys.exit(1)

    c=time.time()
    tag_dict = read_tag(tag,maxlen,minlen)
    print "read tag ",time.time()-c
    c=time.time()
#    tag_dict=makeloci(ps,pe,ns,ne,tagpos)
#    print "combine tag ",time.time()-c
    b=time.time()
    peak_dict = read_interval(peak,span,orient)
    print "read interval ",time.time()-b
    e=time.time()
    #result = distance_length(peak_dict,ps,pe,ns,ne,tagpos,span,maxlen,minlen)
    #print "caculate ",time.time()-e
    result_dict = distance_length(peak_dict,tag_dict,tagpos,span,maxlen,minlen)
    print "caculate ",time.time()-e
    # print result_dict[0]
    output(result_dict,out,minlen,maxlen,span)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

