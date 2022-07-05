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
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader

import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def readBG(bgmatrix):
    pBG = {}
    nBG = {}
    inf = open(bgmatrix)
    for line in inf:
        ll = line.split()
        name = ll[0]
        pBG[name] = float(ll[1])
        nBG[name] = float(ll[2])
    inf.close()
    return pBG,nBG
def sitepro_scan(peak,out_bed,out_matrix_p,out_matrix_n,w_plus,w_minus,bgmatrix,span,gen,lflank,rflank):
    nmer = lflank + rflank
    genome = twobitreader.TwoBitFile(gen)
    pBG,nBG = readBG(bgmatrix)
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    ### dict for distribution
    p_dis ={}
    #p_dis_total = {}
    n_dis = {}
   # n_dis_total = {}
    cleavage_max_p = 0
    cleavage_max_n = 0
    
    outf = open(out_bed,'w')
    for line in inf:### chr start end name motifscore strand FP DNase chip
        ll = line.split()#####  3 below is flanking length
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        ## remove overflow
        if start - span -lflank <= 0:
            continue
        ## get cleavage
        p_sum = list(w_plus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        n_sum = list(w_minus_H.summarize(chrm,start-span,end+span,end-start+2*span).sum_data)
        ## get seqbias
        seq = genome[chrm][(start - span - lflank):(end + span + rflank)]
        if 'N' in seq.upper():
            continue
        pseq = seq[:-1]
        nseq = seq[1:]
        p=[]
        n=[]
        for k in range(len(pseq)  +1 - nmer):
            p.append(pBG[pseq[k:(k+nmer)].upper()])
            n.append(nBG[nseq[k:(k+nmer)].upper()])
        ## get predicted seqbias
        p_assign = []
        n_assign = []
        for bp in range(len(p_sum)- 2*span):
            ptotal = sum(p_sum[bp:(bp+2*span)])
            ntotal = sum(n_sum[bp:(bp+2*span)])
            pbias_per = p[bp+span]*1.0/sum(p[bp:(bp+2*span)])
            nbias_per = n[bp+span]*1.0/sum(n[bp:(bp+2*span)])
            paraw = pbias_per*ptotal
            naraw = nbias_per*ntotal
            p_assign.append(paraw)
            n_assign.append(naraw)
 
            ## calculate real cut distribution
            pc = int(p_sum[bp+span])
            nc = int(n_sum[bp+span])

            if paraw == 0:
                pa = -1
            elif paraw < 1.0/32:
                pa = 0.0
            elif paraw < 1.0/16:
                pa = 1.0/32
            elif paraw < 1.0/8:
                pa = 1.0/16
            elif paraw < 1.0/4:
                pa = 1.0/8
            elif paraw < 1.0/2:
                pa = 1.0/4
            elif paraw < 1.0:
                pa = 1.0/2
            else:
                pa = int(paraw)
            if naraw == 0:
                na = -1
            elif naraw < 1.0/32:
                na = 0.0
            elif naraw < 1.0/16:
                na = 1.0/32
            elif naraw < 1.0/8:
                na = 1.0/16
            elif naraw < 1.0/4:
                na = 1.0/8
            elif naraw < 1.0/2:
                na = 1.0/4
            elif naraw < 1.0:
                na = 1.0/2
            else:
                na = int(naraw)          

            if not p_dis.has_key(pa) :
                p_dis[pa] = {}
            #    p_dis_total[pa] = 0
            if not n_dis.has_key(na) :
                n_dis[na] = {}
            #    n_dis_total[na] = 0
            
            #p_dis_total[pa] += 1
            #n_dis_total[na] += 1
            
            if not p_dis[pa].has_key(pc):
                p_dis[pa][pc] = 0
            if not n_dis[na].has_key(nc):
                n_dis[na][nc] = 0
            p_dis[pa][pc] += 1
            n_dis[na][nc] += 1
        
            cleavage_max_p = max(pc,cleavage_max_p)
            cleavage_max_n = max(nc,cleavage_max_n)
        ### write  real cleavage , seqbias , seqbias predicted cleavage
        newll = ll + p_sum[span:(len(p_sum)-span)] + n_sum[span:(len(n_sum)-span)] + p[span:(len(p)-span)] + n[span:(len(n)-span)] + p_assign + n_assign
        outf.write("\t".join(map(str,newll))+"\n")

    outf.close()
    inf.close()
    ### write real cut distribution matrix
    outfp = open(out_matrix_p,'w')
    outfn = open(out_matrix_n,'w')
    outfp.write("\t".join( ['predict']+['real'+str(i) for i in range(int(cleavage_max_p)+1)] ) + "\n")
    outfn.write("\t".join( ['predict']+['real'+str(i) for i in range(int(cleavage_max_n)+1)] ) + "\n")

    for predict_p in sorted(p_dis.keys()):
        realcuts_p = [0]*(cleavage_max_p+1)
        for p_real in p_dis[predict_p].keys():
            realcuts_p[p_real] += p_dis[predict_p][p_real]
        outfp.write("\t".join(map(str,[predict_p] + realcuts_p))+"\n") 
    outfp.close()

    for predict_n in sorted(n_dis.keys()):
        realcuts_n = [0]*(cleavage_max_n+1)
        for n_real in n_dis[predict_n].keys():
            realcuts_n[n_real] += n_dis[predict_n][n_real]
        outfn.write("\t".join(map(str,[predict_n] + realcuts_n))+"\n") 
    outfn.close()
 
def read_distribution_dict(distribution,pre_uplim,cut_uplim):
    inf = open(distribution)
    dis_dict = {}
    
    for line in inf:
        if line.startswith('predict'):
            continue
        ll = line.split()
        if float(ll[0]) > pre_uplim:
            #print ll[0]
            continue
        dis_dict[float(ll[0])] = map(int,ll[1:(1+cut_uplim)])
    inf.close()
    print 1
    return dis_dict

def pvalue_dict_raw(dis_dict):
    #sdnorm = scipy.stats.distributions.norm_gen()
    p_right = {}
    p_left = {}
    #print dis_dict[-1]
    for pred in sorted(dis_dict.keys()):
        #print pred
        p_right[pred] = []
        p_left[pred] = []
        
        for real in range(len(dis_dict[pred])):
            if sum(dis_dict[pred]) < 10:
                pR = 0.5
                pL = 0.5
            else:
                if real < len(dis_dict[pred]) - 1:
                    pR = (dis_dict[pred][(real)]*1.0+sum(dis_dict[pred][(real+1):])*1.0)/sum(dis_dict[pred])
                else:
                    pR = dis_dict[pred][(real)]*1.0/sum(dis_dict[pred])
                pL = (dis_dict[pred][(real)]*1.0 + sum(dis_dict[pred][:(real)])*1.0)/sum(dis_dict[pred])
            #pval = p# min(p,1-p) *2
            p_right[pred].append(pR)
            p_left[pred].append(pL)
    
            #if pval == 1:
            #    print pval,sdnorm.isf(pval),pred,dis_dict[pred]
            #print pval,z
            #sys.exit(1)
    return p_left,p_right

#def combine_z_to_p(zlist,sdnorm):
#    Z = sum(zlist)*1.0/math.sqrt(len(zlist))
#    P = sdnorm.cdf(Z)
#    Pval = 1-P
#    return Pval
def combine_p_to_chi(plist):
    PI = 0
    for p in plist:
        PI += math.log(p)
    Xsq = -2* PI
    return Xsq
    
def cleavage_pvalue(data,matrix_p,matrix_n,span,out,prelim,cutlim):
    chi = scipy.stats.distributions.chi_gen()
    p_left_plus,p_right_plus = pvalue_dict_raw(read_distribution_dict(matrix_p,prelim,cutlim))
    p_left_minus,p_right_minus = pvalue_dict_raw(read_distribution_dict(matrix_n,prelim,cutlim))
   # for i in range(5):
   #     for j in range(i+1):
   #         print i,j,zdict_p[i][j]
    print matrix_p,matrix_n,'read ok'
    inf = open(data)
    outf = open(out,'w')
    for line in inf:
        ll = line.split()
        each = (len(ll)-5)/6
        p_cut = ll[5:(5+each)]
        n_cut = ll[(5+each):(5+each*2)]
        p_assign = ll[(5+each*4):(5+each*5)]
        n_assign = ll[(5+each*5):(5+each*6)]
        p_left_pval = []
        n_left_pval = []
        p_right_pval = []
        n_right_pval = []
        for i in range(len(p_cut)):
            pc = int(float(p_cut[i]))
            nc = int(float(n_cut[i]))
            paraw = (float(p_assign[i]))
            naraw = (float(n_assign[i]))
            if paraw == 0:
                pa = -1
            elif paraw < 1.0/32:
                pa = 0.0
            elif paraw < 1.0/16:
                pa = 1.0/32
            elif paraw < 1.0/8:
                pa = 1.0/16
            elif paraw < 1.0/4:
                pa = 1.0/8
            elif paraw < 1.0/2:
                pa = 1.0/4
            elif paraw < 1.0:
                pa = 1.0/2
            else:
                pa = int(paraw)
            if naraw == 0:
                na = -1
            elif naraw < 1.0/32:
                na = 0.0
            elif naraw < 1.0/16:
                na = 1.0/32
            elif naraw < 1.0/8:
                na = 1.0/16
            elif naraw < 1.0/4:
                na = 1.0/8
            elif naraw < 1.0/2:
                na = 1.0/4
            elif naraw < 1.0:
                na = 1.0/2
            else:
                na = int(naraw)          

            #print len(zdict_p[pa]),pc
            #sys.exit(1)
            if p_left_plus.has_key(pa) and pc <= len(p_left_plus[pa]) -1 : 
                #print len(zdict_p[pa]),pc
                #sys.exit(1)
                p_left_pval.append(p_left_plus[pa][pc])
                p_right_pval.append(p_right_plus[pa][pc])
                #print pa,pc,zdict_p[pa][pc]
            else:
                p_left_pval.append(1)
                p_right_pval.append(1)
            if p_left_minus.has_key(na) and nc <= len(p_left_minus[na]) -1 : 
                n_left_pval.append(p_left_minus[na][nc])
                n_right_pval.append(p_right_minus[na][nc])
                #print na,nc,zdict_n[na][nc]
            else:
                n_left_pval.append(1)
                n_right_pval.append(1)
            #n_zscore.append(zdict_n[na][nc])
        Xs_left = []
        Ps_left = []
        Xs_right = []
        Ps_right = []
        #Xs_LR = []
        #Ps_LR = []
        #Xs_RL = []
        #Ps_RL = []
        for i in range(len(p_left_pval) - 2*span):
            X_left = combine_p_to_chi(p_left_pval[i:(i+2*span)] + n_left_pval[i:(i+2*span)] )
            P_left = 1-chi.cdf(math.sqrt(X_left),8*span)
            #P_n = combine_z_to_p(n_zscore[i:(i+2*span)],sdnorm)
            Ps_left.append(P_left)
            Xs_left.append(X_left)
            X_right = combine_p_to_chi(p_right_pval[i:(i+2*span)] + n_right_pval[i:(i+2*span)] )
            P_right = 1-chi.cdf(math.sqrt(X_right),8*span)
            Ps_right.append(P_right)
            Xs_right.append(X_right)
            
            #x_lr = combine_p_to_chi(p_left_pval[i:(i+span)] + n_left_pval[i:(i+span)] + p_right_pval[(i+span):(i+span*2)] + n_right_pval[(i+span):(i+span*2)])
            #x_rl = combine_p_to_chi(p_right_pval[i:(i+span)] + n_right_pval[i:(i+span)] + p_left_pval[(i+span):(i+span*2)] + n_left_pval[(i+span):(i+span*2)])
            #p_lr = 1-chi.cdf(math.sqrt(x_lr),8*span)
            #p_rl = 1-chi.cdf(math.sqrt(x_rl),8*span)
            #Xs_LR.append(x_lr)
            #Xs_RL.append(x_rl)
            #Ps_LR.append(p_lr)
            #Ps_RL.append(p_rl)
        newll = ll + p_left_pval + n_left_pval + p_right_pval + n_right_pval  + Ps_left  + Ps_right #+  Ps_LR + Ps_RL #+ n_Pval
        outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
    inf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="if raw , -i is 5 column summit200bp file\
                               if predict, -i is 5 column summit200bp + pcut,ncut,pbias,nbias,passign,nassign, 400 column each")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("--w1",dest="w_plus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_plus_50_100.bw",
                         help="")
    optparser.add_option("--w2",dest="w_minus",type="str",default = "/home/sh430/Project/ChongzhiDNase/Data/bwTrack/Lncap_DHT/Lncap_DNase_minus_50_100.bw",
                         help="")
    optparser.add_option("-b","--bgmatrix",dest="bgmatrix",type="str",
                         help="sequence bias matrix")
    optparser.add_option("-m","--mode",dest="mode",type="str",default = "raw",
                         help="mode , choose from raw and predict, \
                                      if raw, input -i -o --w1/2 -b --genome\
                                      if predict , input -i -o --m1/2 ")
    optparser.add_option("--m1",dest="m_plus",type="str",
                         help="used if predict mode , 1 is plus , 2 is minus")
    optparser.add_option("--m2",dest="m_minus",type="str",
                         help="")                                                                          
#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 10,
                         help="region for get total signal in single bp, default = 10 means +-10bp(total 20bp) signal as total for each bp")
    optparser.add_option("--Zspan",dest="Zspan",type="int",default = 10,
                         help="region for combine Zscore, default = 10 means +-10bp(total 20bp) signal as total for each bp")

    optparser.add_option("--genome",dest="genome",type="str",
                         help="2bit format")
    optparser.add_option("--left",dest="leftflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--right",dest="rightflank",type="int",default = 3,
                         help="flnaking region for seqbias , 6-mer means left=right=3")
    optparser.add_option("--prelim",dest="prelim",type="int",default = 500,
                         help="preidct cut up limit , any predicted number greater than this will get Zscore as 0")
    optparser.add_option("--cutlim",dest="cutlim",type="int",default = 5000,
                         help="real cut up limit , any predicted number greater than this will get Zscore as 0")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)


    interval = options.interval
    out_bed = options.output + '_predict.bed'
    out_matrix_p = options.output + '_Pdistribution.txt'
    out_matrix_n = options.output + '_Ndistribution.txt'
    out = options.output + '_pval.bed'
    w_plus = options.w_plus
    w_minus = options.w_minus
    m_plus = options.m_plus
    m_minus = options.m_minus
    bgmatrix = options.bgmatrix
    gen = options.genome
    lflank = options.leftflank
    rflank = options.rightflank
    
    if options.mode == 'raw':
        sitepro_scan(interval,out_bed,out_matrix_p,out_matrix_n,w_plus,w_minus,bgmatrix,options.Cspan,gen,lflank,rflank)
        cleavage_pvalue(out_bed,out_matrix_p,out_matrix_n,options.Zspan,out,options.prelim,options.cutlim)
    elif options.mode == 'predict':
        out_bed = interval
        out_matrix_p = m_plus
        out_matrix_n = m_minus
        cleavage_pvalue(out_bed,out_matrix_p,out_matrix_n,options.Zspan,out,options.prelim,options.cutlim)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

