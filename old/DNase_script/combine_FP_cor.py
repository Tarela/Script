#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:
for new motif+TC+FOS+cleavage+seqbias output file , collect FOS (column 8) and correlation coefficient for cleavage ~ seqbias(central 50bp)

"""
import os,sys
import scipy.stats.stats as sat
import numpy
def analysisFile(inputfile,outf):
    FOSlist = []
    inf = open(inputfile)
    testline = inf.readline().split()[8:]
#    print len(testline)
    ml = (len(testline))/3
    all_c = numpy.array([0]*ml)
    all_b = numpy.array([0]*ml)
    count=0
    name = inputfile.lstrip('.bed')
    
    for line in inf:
        count+=1
        ll = line.split()
        FOSlist.append(float(ll[7]))
        bpnumber = (len(ll)-8)/6
#        print bpnumber
        #print f
        #print len(ll)
        pc = ll[8:(8+bpnumber)]
        nc = ll[(8+bpnumber):(8+bpnumber*2)]
        pI = ll[(8+bpnumber*2):(8+bpnumber*3)]
        nI = ll[(8+bpnumber*3):(8+bpnumber*4)]
        pb = ll[(8+bpnumber*4):(8+bpnumber*5)]
        nb = ll[(8+bpnumber*5):(8+bpnumber*6)]
#        print len(pc)
        c = numpy.array(map(float,pc+nc))
        b = numpy.array(map(float,pb+nb))
#        print len(c)
        #cc = sat.pearsonr(c,b)[0]
        #print c
        #print b
        all_c += c
        all_b += b
    
    ave_c = all_c*1.0/count
    ave_b = all_b*1.0/count
    
    c_c = sat.pearsonr(ave_c,ave_b)[0]
    newll = [name,c_c, (numpy.mean(FOSlist))]# + [-1]*(5000-count)
 #   print len(FOSlist)
 #   print 5000-count
    outf.write("\t".join(map(str,newll))+"\n")
    #outf.close()
    inf.close()
        
allfile = os.listdir(sys.argv[1])
allfile.sort()
outf = open(sys.argv[2],'w')
for f in allfile:
    if f.endswith(sys.argv[3]):
        analysisFile(f,outf)
        print f,'done'
outf.close()
