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

def makeoverlap(curlist):
    a=[0]*5
    for i in curlist:
        for j in range(1,6):
            if int(i[j]) > 0 :
                a[j-1] += 1
    return a

def single_binplot(inputlist,binsize):
    summary = [[],[],[],[],[]]
    for binnum in range(len(inputlist)/binsize):
        current_list = inputlist[binnum*binsize:(binnum+1)*binsize]
        overlap = makeoverlap(current_list)
        for i in range(5):
            summary[i].append(overlap[i])
    return summary

def binplot(input_table,out,binsize,sortby):
    DHSlist = []
    motiflist = []
    nonmotiflist = []
    inf = open(input_table)
    for line in inf:
        ll = line.split()
        newll = [float(ll[sortby])] + ll[17:22]
        DHSlist.append(newll)
        if int(ll[22] ) > 0 :
            motiflist.append(newll)
        else:
            nonmotiflist.append(newll)
    inf.close()
    DHSlist.sort( key=lambda x:x[0],reverse=True )
    motiflist.sort( key=lambda x:x[0],reverse=True )
    nonmotiflist.sort( key=lambda x:x[0],reverse=True )

    DHSplot = single_binplot(DHSlist,binsize)
    motifplot = single_binplot(motiflist,binsize)
    nonmotifplot = single_binplot(nonmotiflist,binsize)

    outf = open(out+".r",'w')
    outf.write('''
oldAR_all<-c(%s)
C19_all<-c(%s)
H280_all<-c(%s)
N20_all<-c(%s)
combine_all<-c(%s)
oldAR_yes<-c(%s)
C19_yes<-c(%s)
H280_yes<-c(%s)
N20_yes<-c(%s)
combine_yes<-c(%s)
oldAR_no<-c(%s)
C19_no<-c(%s)
H280_no<-c(%s)
N20_no<-c(%s)
combine_no<-c(%s)
pdf(file="%s")
par(mfrow=c(3,1),mar=c(2,4,4,2))
plot(combine_all,type="l",ylim=c(0,400),ylab="number in 500 overlap ARpeak",main="DHS overlap with all AR, rank high -> low")
lines(C19_all,col="red")
lines(H280_all,col="blue")
lines(N20_all,col="green")
lines(oldAR_all,col="orange")

plot(combine_yes,type="l",ylim=c(0,400),ylab="number in 500 overlap ARpeak",main="DHS overlap with AR w/ motif")
lines(C19_yes,col="red")
lines(H280_yes,col="blue")
lines(N20_yes,col="green")
lines(oldAR_yes,col="orange")

plot(combine_no,type="l",ylim=c(0,400),ylab="number in 500 overlap ARpeak",main="DHS overlap with AR w/o motif")
lines(C19_no,col="red")
lines(H280_no,col="blue")
lines(N20_no,col="green")
lines(oldAR_no,col="orange")
dev.off()
'''%(str(DHSplot[0])[1:-1],str(DHSplot[1])[1:-1],str(DHSplot[2])[1:-1],str(DHSplot[3])[1:-1],str(DHSplot[4])[1:-1],str(motifplot[0])[1:-1],str(motifplot[1])[1:-1],str(motifplot[2])[1:-1],str(motifplot[3])[1:-1],str(motifplot[4])[1:-1],str(nonmotifplot[0])[1:-1],str(nonmotifplot[1])[1:-1],str(nonmotifplot[2])[1:-1],str(nonmotifplot[3])[1:-1],str(nonmotifplot[4])[1:-1],out+".pdf")) 
    outf.close()
    os.system('Rscript %s'%(out+".r"))
    #print sum(DHSplot[0]),sum(DHSplot[1]),sum(DHSplot[2]),sum(DHSplot[3]),sum(DHSplot[4])
    #for i in range(5):
    #    print DHSplot[i]

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input_table",dest="input_table",type="str",default = "Lncap_DHS_ARpeak_motif.bed",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="")
    optparser.add_option("-b","--binsize",dest="binsize",type="int",default=500,
                         help="")
    optparser.add_option("-s","--sortby",dest="sortby",type="int",default=13,
                         help="")

#========minor options=============

    (options,args) = optparser.parse_args()

    input_table = options.input_table
    out = options.output
    binsize = options.binsize
    sortby  = options.sortby
    if not input_table:
        optparser.print_help()
        sys.exit(1)
    
    binplot(input_table,out,binsize,sortby)

    


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

